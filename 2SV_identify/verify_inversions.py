#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_inversions.py

Efficient inversion verification for master_table:
 - split-read (SA tag) and soft-clip support per breakpoint
 - mean coverage inside inversion vs flanks (depth_ratio)
 - orientation consistency inside inversion
 - faster by avoiding per-base loops
 - merges PASS + WEAK into PASS for simplicity

"""

import argparse
import pandas as pd
import pysam
import sys


def parse_inv_row(row):
    """Parse inversion row, supporting both named and unnamed columns."""
    if all(c in row.index for c in ['chrom', 'start', 'end']):
        chrom = str(row['chrom'])
        start = int(row['start'])
        end = int(row['end'])
        length = int(row['length']) if 'length' in row.index else (end - start + 1)
    else:
        chrom = str(row.iloc[0])
        start = int(row.iloc[1])
        end = int(row.iloc[2])
        length = int(row.iloc[3]) if len(row) > 3 else (end - start + 1)
    return chrom, start, end, length


def count_breakpoint_support(bam, chrom, pos1, window=100, min_mapq=20, min_clip=20):
    """Count split-read (SA) and soft-clipped reads around a breakpoint."""
    split_count, soft_count, total_reads = 0, 0, 0
    start0 = max(0, pos1 - window - 1)
    end0 = pos1 + window
    try:
        for read in bam.fetch(chrom, start0, end0):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            total_reads += 1
            if read.has_tag("SA"):
                sa_entries = [x for x in read.get_tag("SA").split(";") if x]
                for entry in sa_entries:
                    fields = entry.split(",")
                    if len(fields) >= 3 and fields[0] == chrom:
                        split_count += 1
                        break
            if read.cigartuples:
                if any(op == 4 and length >= min_clip for op, length in read.cigartuples):
                    soft_count += 1
    except ValueError:
        return {'split': 0, 'soft': 0, 'total': 0}
    return {'split': split_count, 'soft': soft_count, 'total': total_reads}


def mean_coverage_region(bam, chrom, start1, end1):
    """Compute mean coverage across region [start1,end1]."""
    if end1 < start1:
        return 0.0
    start0, end0 = start1 - 1, end1
    try:
        cov = bam.count_coverage(chrom, start=start0, end=end0, quality_threshold=0)
    except ValueError:
        return 0.0
    if not cov or len(cov[0]) == 0:
        return 0.0
    depth_array = [sum(x) for x in zip(*cov)]
    return float(sum(depth_array) / len(depth_array))


def orientation_consistency(bam, chrom, start1, end1, min_mapq=20, max_reads=50000):
    """Proportion of reads in majority orientation inside inversion."""
    fwd, rev, n = 0, 0, 0
    try:
        for read in bam.fetch(chrom, start1 - 1, end1):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            if read.is_reverse:
                rev += 1
            else:
                fwd += 1
            n += 1
            if n > max_reads:  # cap to avoid huge inversions being slow
                break
    except ValueError:
        return 0.0
    total = fwd + rev
    return max(fwd, rev) / total if total > 0 else 0.0


def main():
    parser = argparse.ArgumentParser(description="Verify inversions from BAM and produce support metrics for master_table")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (indexed)")
    parser.add_argument("-i", "--inv", required=True, help="Inversion CSV/TSV (chrom,start,end,length)")
    parser.add_argument("-o", "--out", required=True, help="Output CSV/TSV file")
    parser.add_argument("--window", type=int, default=100, help="Window around breakpoints [100]")
    parser.add_argument("--flank", type=int, default=1000, help="Flank size for coverage comparison [1000]")
    parser.add_argument("--min_mapq", type=int, default=20, help="Min MAPQ to count reads [20]")
    parser.add_argument("--min_clip", type=int, default=20, help="Min soft-clip length [20]")
    parser.add_argument("--split_thresh", type=int, default=1, help="Min split reads per side [1]")
    parser.add_argument("--orient_thresh", type=float, default=0.5, help="Min orientation consistency [0.5]")
    parser.add_argument("--sep", default=",", help="Separator of input CSV [comma]")
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    inv_df = pd.read_csv(args.inv, sep=args.sep)
    results = []

    for _, row in inv_df.iterrows():
        chrom, start, end, length = parse_inv_row(row)
        if chrom not in bam.references:
            results.append([chrom, start, end, length, 0, 0, 0, 0, 0, 0, 0.0, 0.0, None, 0.0, "FAIL"])
            continue
        left = count_breakpoint_support(bam, chrom, start, args.window, args.min_mapq, args.min_clip)
        right = count_breakpoint_support(bam, chrom, end, args.window, args.min_mapq, args.min_clip)
        cov_inside = mean_coverage_region(bam, chrom, start, end)
        cov_left = mean_coverage_region(bam, chrom, max(1, start - args.flank), start - 1)
        cov_right = mean_coverage_region(bam, chrom, end + 1, min(bam.get_reference_length(chrom), end + args.flank))
        flank_mean = (cov_left + cov_right) / 2 if (cov_left > 0 and cov_right > 0) else max(cov_left, cov_right)
        depth_ratio = cov_inside / flank_mean if flank_mean > 0 else None
        orient = orientation_consistency(bam, chrom, start, end, args.min_mapq)

        # unified decision: PASS if both breakpoints have >= split_thresh and orient >= orient_thresh
        status = "PASS" if (left['split'] >= args.split_thresh and right['split'] >= args.split_thresh and orient >= args.orient_thresh) else "FAIL"

        results.append([chrom, start, end, length,
                        left['split'], right['split'],
                        left['soft'], right['soft'],
                        left['total'], right['total'],
                        cov_inside, flank_mean, depth_ratio,
                        orient, status])

    out_df = pd.DataFrame(results, columns=[
        "chrom", "start", "end", "length",
        "left_split", "right_split",
        "left_soft", "right_soft",
        "left_total", "right_total",
        "coverage_inside_mean", "coverage_flank_mean", "depth_ratio",
        "CONSIST_ORIENT", "status"
    ])
    out_df.to_csv(args.out, index=False)
    bam.close()
    print(f"[done] Results saved to {args.out}")


if __name__ == "__main__":
    main()