#!/usr/bin/env python3
"""
summarize_inversions_genes.py

  --inversions All_PASS_inversions.bed  # bed-like: chrom start end sampleID (whitespace separated)
  --gtf-dir /path/to/gtf_folder         # folder containing {sample}.braker.gtf (or use --gtf-pattern)
  --downstream 1000     
  --gtf-pattern "{sample}.braker.gtf"
  --out-prefix "inv_genes"

  {out_prefix}_per_inversion.tsv
  {out_prefix}_per_sample.tsv
"""
import os, sys, argparse, re, glob, bisect, csv
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Summarize genes inside / at breakpoints / downstream for inversions.")
    p.add_argument("--inversions", required=True, help="All_PASS_inversions.bed (chrom start end sampleID)")
    p.add_argument("--gtf-dir", required=True, help="Directory with GTF files")
    p.add_argument("--downstream", type=int, default=1000, help="Downstream window in bp (default 1000)")
    p.add_argument("--gtf-pattern", default="{sample}.braker.gtf", help="GTF filename pattern containing {sample}")
    p.add_argument("--out-prefix", default="inv_genes", help="Output prefix")
    return p.parse_args()

# robust attribute parser for GTF 9th column
attr_kv_re = re.compile(r'(\S+)\s+"([^"]+)"')
def parse_attributes(attr_str):
    d = {}
    for part in attr_str.strip().split(';'):
        part = part.strip()
        if not part:
            continue
        m = attr_kv_re.match(part)
        if m:
            d[m.group(1)] = m.group(2)
        else:
            if '=' in part:
                k,v = part.split('=',1)
                d[k.strip()] = v.strip()
            else:
                # fallback
                d[part] = part
    return d

def parse_gtf(gtf_path):
    """
    Parse GTF and return index:
      index[chrom] = {
         'intervals': [(start,end,gene_id), ...] sorted by start,
         'starts': [start,...],
         'ends': [end,...]
      }
    For each gene_id we take min(start) and max(end) across all lines with same gene_id on same chrom.
    """
    gene_coords = {}  # key = (chrom,gene_id) -> [minstart,maxend]
    try:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith('#'): continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9: continue
                chrom = cols[0]
                try:
                    start = int(cols[3])
                    end = int(cols[4])
                except:
                    continue
                attrs = parse_attributes(cols[8])
                gene_id = attrs.get('gene_id') or attrs.get('gene') or attrs.get('ID') or attrs.get('transcript_id')
                if not gene_id:
                    continue
                key = (chrom, gene_id)
                if key in gene_coords:
                    gene_coords[key][0] = min(gene_coords[key][0], start)
                    gene_coords[key][1] = max(gene_coords[key][1], end)
                else:
                    gene_coords[key] = [start, end]
    except Exception as e:
        raise RuntimeError(f"Failed to parse GTF {gtf_path}: {e}")

    index = {}
    for (chrom,gid),(s,e) in gene_coords.items():
        index.setdefault(chrom, []).append((s,e,gid))
    # sort and create helper lists
    for chrom in list(index.keys()):
        lst = sorted(index[chrom], key=lambda x: x[0])
        starts = [x[0] for x in lst]
        ends = [x[1] for x in lst]
        index[chrom] = {'intervals': lst, 'starts': starts, 'ends': ends}
    return index

def find_overlaps(index, chrom, a, b):
    """Return list of gene_ids overlapping [a,b] (inclusive).
       If chrom not in index -> empty list.
    """
    if chrom not in index:
        return []
    starts = index[chrom]['starts']
    intervals = index[chrom]['intervals']
    # all intervals with start <= b are in intervals[:i]
    i = bisect.bisect_right(starts, b)
    res = []
    for j in range(i):
        s,e,gid = intervals[j]
        if e >= a:
            res.append(gid)
    # remove duplicates but keep order
    return list(dict.fromkeys(res))

def find_genes_at_pos(index, chrom, pos):
    return find_overlaps(index, chrom, pos, pos)

def find_gtf_for_sample(gtf_dir, pattern, sample):
    # first try direct formatting
    candidate = os.path.join(gtf_dir, pattern.format(sample=sample))
    if os.path.exists(candidate):
        return candidate
    # fallback: any file in gtf_dir containing sample string and ending with .gtf or .gtf.gz or .braker.gtf
    for ext in ("*.gtf","*.gtf.gz","*.braker*","*.gtf*"):
        matches = glob.glob(os.path.join(gtf_dir, f"*{sample}*"))
        if matches:
            # prefer .gtf if exact match
            for m in matches:
                if m.endswith('.gtf'):
                    return m
            return matches[0]
    return None

def main():
    args = parse_args()
    inv_file = args.inversions
    gtf_dir = args.gtf_dir
    window = args.downstream
    gtf_pattern = args.gtf_pattern
    out_prefix = args.out_prefix

    # read inversions -> group by sample
    sample_to_invs = defaultdict(list)
    with open(inv_file) as fh:
        for line in fh:
            if line.strip() == "" or line.startswith('#'): continue
            cols = line.strip().split()
            if len(cols) < 3: continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except:
                continue
            sample = cols[3] if len(cols) >= 4 else "unknown"
            sample_to_invs[sample].append((chrom, start, end))

    out_inv = out_prefix + "_per_inversion.tsv"
    out_sample = out_prefix + "_per_sample.tsv"

    # prepare output header
    header = [
        "sample","chrom","inv_start","inv_end","inv_id",
        "genes_inside_count","genes_inside_ids",
        "left_bp_genes_count","left_bp_genes_ids",
        "right_bp_genes_count","right_bp_genes_ids",
        f"downstream_{window}bp_genes_count", f"downstream_{window}bp_genes_ids"
    ]

    # per-sample accumulators for summary
    sample_summaries = {}

    with open(out_inv, 'w', newline='') as outf:
        w = csv.writer(outf, delimiter='\t')
        w.writerow(header)

        for sample in sorted(sample_to_invs.keys()):
            invs = sample_to_invs[sample]
            # find gtf
            gtf_path = find_gtf_for_sample(gtf_dir, gtf_pattern, sample)
            if gtf_path:
                try:
                    index = parse_gtf(gtf_path)
                except Exception as e:
                    sys.stderr.write(f"[WARN] Failed parsing GTF for {sample}: {e}\n")
                    index = {}
            else:
                sys.stderr.write(f"[WARN] No GTF found for sample {sample} (tried pattern {gtf_pattern}). All results will be 0.\n")
                index = {}

            # summary accumulators
            unique_inside = set()
            unique_bp = set()
            unique_down = set()
            inversions_with_inside = 0
            total_inversions = len(invs)
            total_genes_inside_count = 0  # counting duplicates across inversions
            for (chrom,start,end) in invs:
                inv_id = f"{chrom}:{start}-{end}"
                inside = find_overlaps(index, chrom, start, end)
                left_bp = find_genes_at_pos(index, chrom, start)
                right_bp = find_genes_at_pos(index, chrom, end)
                down_start = end + 1
                down_end = end + window
                downstream = find_overlaps(index, chrom, down_start, down_end) if window > 0 else []

                if inside:
                    inversions_with_inside += 1
                total_genes_inside_count += len(inside)
                for g in inside: unique_inside.add(g)
                for g in (left_bp + right_bp): unique_bp.add(g)
                for g in downstream: unique_down.add(g)

                def fmt(lst):
                    return ";".join(lst) if lst else "0"

                row = [
                    sample, chrom, start, end, inv_id,
                    len(inside), fmt(inside),
                    len(left_bp), fmt(left_bp),
                    len(right_bp), fmt(right_bp),
                    len(downstream), fmt(downstream)
                ]
                w.writerow(row)

            sample_summaries[sample] = {
                'total_inversions': total_inversions,
                'inversions_with_genes_inside': inversions_with_inside,
                'total_genes_inside_count': total_genes_inside_count,
                'unique_genes_inside': len(unique_inside),
                'unique_bp_genes': len(unique_bp),
                'unique_down_genes': len(unique_down)
            }

    # write per-sample summary
    with open(out_sample, 'w', newline='') as sh:
        sw = csv.writer(sh, delimiter='\t')
        sh_header = [
            "sample","total_inversions","inversions_with_genes_inside",
            "total_genes_inside_count","unique_genes_inside",
            "unique_breakpoint_genes","unique_downstream_genes"
        ]
        sw.writerow(sh_header)
        for sample, data in sorted(sample_summaries.items()):
            sw.writerow([
                sample,
                data['total_inversions'],
                data['inversions_with_genes_inside'],
                data['total_genes_inside_count'],
                data['unique_genes_inside'],
                data['unique_bp_genes'],
                data['unique_down_genes']
            ])
    print(f"[DONE] per-inversion table -> {out_inv}")
    print(f"[DONE] per-sample summary -> {out_sample}")

if __name__ == "__main__":
    main()