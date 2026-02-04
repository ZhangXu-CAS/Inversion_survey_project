#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_inversion_pi.py
For each sample:
 - Compute π for inversion regions (from BED)
 - Compute π for scaffold non-inversion regions
 - Compute ratio = π_inv / π_noninv
"""

import os, sys, argparse, pysam, pandas as pd, numpy as np

def safe_makedirs(path):
    if not os.path.exists(path):
        os.makedirs(path)

def compute_divergence_from_bam(bam, chrom, start0, end1):
    """Compute pairwise π from BAM for region"""
    try:
        iter_reads = bam.fetch(chrom, start0, end1)
    except ValueError:
        return None

    aligned_ref_len = mismatches = 0
    for read in iter_reads:
        if read.is_unmapped or read.is_secondary: continue
        try:
            nm = read.get_tag('NM')  # total mismatches + indels
        except KeyError:
            nm = 0
        aligned_ref_len += read.query_alignment_length
        mismatches += nm
    if aligned_ref_len == 0:
        return None
    return mismatches / aligned_ref_len

def merge_intervals(df):
    """Merge overlapping intervals in BED DataFrame"""
    merged = []
    for chrom, group in df.groupby('chrom'):
        group = group.sort_values('start0')
        cur_s, cur_e = None, None
        for s, e in zip(group.start0, group.end1):
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((chrom, cur_s, cur_e))
    return merged

def compute_noninv_regions(chrom, chrom_len, inv_intervals):
    """Subtract inversion regions from full scaffold"""
    inv_intervals = sorted(inv_intervals, key=lambda x: x[0])
    noninv = []
    prev_end = 0
    for s, e in inv_intervals:
        if s > prev_end:
            noninv.append((prev_end, s))
        prev_end = max(prev_end, e)
    if prev_end < chrom_len:
        noninv.append((prev_end, chrom_len))
    return noninv

def main(args):
    sample = args.sample
    print(f"[INFO] Computing inversion and non-inversion π for {sample}")

    outdir = args.outdir
    safe_makedirs(outdir)

    bed_df = pd.read_csv(args.bed, sep='\t', header=None, names=['chrom', 'start0', 'end1'])
    bam = pysam.AlignmentFile(args.bam, "rb")

    results = []

    # ---- Compute inversion π ----
    for i, row in bed_df.iterrows():
        chrom, s, e = row.chrom, int(row.start0), int(row.end1)
        pi_inv = compute_divergence_from_bam(bam, chrom, s, e)
        if pi_inv is None:
            print(f"[WARN] No valid alignment for {chrom}:{s}-{e}, skipping")
            continue
        if pi_inv > 0.1:
            print(f"[WARN] Unusually high π={pi_inv:.3f} for {chrom}:{s}-{e}, skipping")
            continue
        results.append(dict(sample=sample, chrom=chrom, start=s, end=e, pi_inv=pi_inv))

    # ---- Compute scaffold-level π (excluding inversions) ----
    chroms = sorted(set(bed_df.chrom))
    scaffold_noninv = []

    for chrom in chroms:
        chrom_len = bam.get_reference_length(chrom)
        inv_intervals = [(int(s), int(e)) for s, e in bed_df.query("chrom == @chrom")[['start0','end1']].values]
        merged = merge_intervals(pd.DataFrame(inv_intervals, columns=['start0','end1']).assign(chrom=chrom))
        merged_intervals = [(s,e) for c,s,e in merged if c==chrom]
        noninv_regions = compute_noninv_regions(chrom, chrom_len, merged_intervals)
        pis = []
        for s, e in noninv_regions:
            pi_val = compute_divergence_from_bam(bam, chrom, s, e)
            if pi_val is not None and pi_val < 0.1:
                pis.append(pi_val)
        if pis:
            pi_noninv = np.mean(pis)
            scaffold_noninv.append((chrom, pi_noninv))

    scaffold_noninv = dict(scaffold_noninv)

    # ---- Combine ----
    for r in results:
        chrom = r['chrom']
        r['pi_noninv'] = scaffold_noninv.get(chrom, np.nan)
        if not np.isnan(r['pi_noninv']):
            r['pi_ratio'] = r['pi_inv'] / r['pi_noninv']
        else:
            r['pi_ratio'] = np.nan

    df = pd.DataFrame(results)
    out_csv = os.path.join(outdir, f"{sample}.inversion_pi.csv")
    df.to_csv(out_csv, index=False)
    print(f"[DONE] Output written: {out_csv}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Compute inversion vs non-inversion divergence (π)")
    ap.add_argument("--bed", required=True, help="BED file of inversions (chrom start end)")
    ap.add_argument("--bam", required=True, help="hap2_vs_hap1 BAM (sorted, indexed)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", default="pi_results")
    args = ap.parse_args()
    main(args)