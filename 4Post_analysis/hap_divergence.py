#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hap_divergence.py

Compute haplotype-pair divergence metrics from minimap2 alignment BAM (hap2 -> hap1).
Outputs:
 - <outprefix>.alignments.csv   : per-alignment stats
 - <outprefix>.windows.tsv      : per-bin stats (chrom, start, end, covered, mismatches, ins, del, divergence, indel_rate)
 - <outprefix>.genome_summary.csv : genome-wide summary (aligned_bases, mismatches, ins, del, pi)
 - optional: <outprefix>.genes.tsv if --gff provided (aggregated by gene)

Requirements: pysam, pandas, numpy
Install: pip install pysam pandas numpy

Usage example:
  # convert SAM to sorted BAM first:
  samtools view -bS aln.sam | samtools sort -o aln.sorted.bam
  samtools index aln.sorted.bam

  python3 hap_divergence.py -b aln.sorted.bam -o sample1.hapdiv --bin 10000 \
      --min_mapq 20 --gff genes.gff --repeat repeatmasker.bed

Notes:
 - Script treats input BAM as assembly-vs-assembly alignments (query contigs -> reference).
 - It's robust to eqx (=,X) CIGAR, or M+MD+NM tags fallback.
"""

import argparse
import re
import math
import sys
from collections import defaultdict
import pysam
import numpy as np
import pandas as pd

# ---------------------- helpers ----------------------
def parse_md_mismatches(mdstr):
    """
    Parse MD string to count mismatches and deletions.
    MD example: '10A5^CT6G2'  -> mismatches letters: A G ; deletions after ^: CT (len2)
    Returns (n_mismatches, n_del_bases)
    """
    if mdstr is None:
        return 0, 0
    # count deletions: occurrences of ^ followed by letters
    del_matches = re.findall(r'\^([A-Z]+)', mdstr)
    del_bases = sum(len(x) for x in del_matches) if del_matches else 0
    # remove deletion parts for mismatches
    md_no_del = re.sub(r'\^[A-Z]+', '', mdstr)
    # mismatches are single letters in between numbers
    mismatches = re.findall(r'[A-Z]', md_no_del)
    return len(mismatches), del_bases

def add_count_to_bins(chrom_bins, chrom_len, bin_size, seg_start0, seg_end0, count):
    """
    Distribute 'count' (an integer or float) that corresponds to reference interval
    [seg_start0, seg_end0) across overlapping bins proportionally to overlap length.
    chrom_bins: numpy array to add to
    """
    if seg_end0 <= seg_start0:
        return
    start_bin = seg_start0 // bin_size
    end_bin = (seg_end0 - 1) // bin_size
    total_len = seg_end0 - seg_start0
    for b in range(start_bin, end_bin + 1):
        bin_start = b * bin_size
        bin_end = min(chrom_len, (b + 1) * bin_size)
        ov = min(seg_end0, bin_end) - max(seg_start0, bin_start)
        if ov > 0:
            chrom_bins[b] += count * (ov / total_len)

def safe_get_tag(read, tag):
    try:
        return read.get_tag(tag)
    except (KeyError, AttributeError):
        return None

# ---------------------- main processing ----------------------
def process_bam(bam_path, outprefix, bin_size=10000, min_mapq=20, max_align_records=None,
                gff_path=None, repeat_bed=None):
    bam = pysam.AlignmentFile(bam_path, "rb")
    # build chrom -> length
    refs = list(bam.references)
    rlen = dict(zip(refs, bam.lengths))
    # initialize bins dict: per chrom arrays for covered, mismatches, ins, del
    chrom_bins = {}
    n_bins = {}
    for chrom in refs:
        nb = (rlen[chrom] + bin_size - 1) // bin_size
        n_bins[chrom] = nb
        chrom_bins[chrom] = {
            'covered': np.zeros(nb, dtype=float),
            'mismatch': np.zeros(nb, dtype=float),
            'ins': np.zeros(nb, dtype=float),
            'del': np.zeros(nb, dtype=float),
            'match': np.zeros(nb, dtype=float)
        }

    # per-alignment rows
    align_rows = []

    total_aligned = 0.0
    total_mism = 0.0
    total_ins = 0.0
    total_del = 0.0

    # iterate over alignments
    nproc = 0
    for read in bam.fetch(until_eof=True):
        nproc += 1
        if max_align_records and nproc > max_align_records:
            break
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        # skip secondary? keep supplementary to capture multi-part alignments
        if read.is_secondary:
            continue

        qname = read.query_name
        rname = read.reference_name
        rstart0 = read.reference_start  # 0-based
        rend0 = read.reference_end      # exclusive
        # aligned_ref_len: count of reference-consuming ops (M, =, X, D)
        # pysam.cigartuples: list of (op,length) where op numeric: see SAM spec
        # codes: 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7/=,8=X
        cig = read.cigartuples
        if cig is None:
            continue
        ref_consumed = 0
        ins_bases = 0
        del_bases = 0
        matches = 0
        mismatches = 0

        # if CIGAR has 7/8 (eqx), prefer using them
        has_eqx = any(op in (7,8) for op, l in cig)

        if has_eqx:
            # handle op7 and op8 explicitly, D consumes reference but isn't match/mismatch
            pos_ref = rstart0
            for op, l in cig:
                if op == 7:  # =
                    # add matches
                    add_count_to_bins(chrom_bins[rname]['match'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    add_count_to_bins(chrom_bins[rname]['covered'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    matches += l
                    ref_consumed += l
                    pos_ref += l
                elif op == 8:  # X mismatch
                    add_count_to_bins(chrom_bins[rname]['mismatch'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    add_count_to_bins(chrom_bins[rname]['covered'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    mismatches += l
                    ref_consumed += l
                    pos_ref += l
                elif op == 2:  # D
                    add_count_to_bins(chrom_bins[rname]['del'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    add_count_to_bins(chrom_bins[rname]['covered'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    del_bases += l
                    ref_consumed += l
                    pos_ref += l
                elif op == 1:  # I - insertion relative to reference: assign at pos_ref
                    # attribute insertion to the bin containing pos_ref (or previous if at start)
                    bin_idx = min((pos_ref)//bin_size, n_bins[rname]-1)
                    chrom_bins[rname]['ins'][bin_idx] += l
                    ins_bases += l
                elif op in (0,):  # M (fallback) treat as covered
                    add_count_to_bins(chrom_bins[rname]['covered'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    ref_consumed += l
                    pos_ref += l
                else:
                    # ignore S/H/P/N
                    if op in (3,): # N
                        ref_consumed += l
                        pos_ref += l
            # compute alignment-level percent identity
        else:
            # no eqx: rely on MD/NM + CIGAR I/D
            # compute ins/del from CIGAR
            pos_ref = rstart0
            for op, l in cig:
                if op == 2:  # D
                    add_count_to_bins(chrom_bins[rname]['del'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    add_count_to_bins(chrom_bins[rname]['covered'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    del_bases += l
                    ref_consumed += l
                    pos_ref += l
                elif op == 1:  # I
                    bin_idx = min((pos_ref)//bin_size, n_bins[rname]-1)
                    chrom_bins[rname]['ins'][bin_idx] += l
                    ins_bases += l
                elif op in (0,7,8):  # M or eqx treat as covered for now
                    add_count_to_bins(chrom_bins[rname]['covered'], rlen[rname], bin_size, pos_ref, pos_ref + l, l)
                    ref_consumed += l
                    pos_ref += l
                else:
                    if op in (3,):
                        ref_consumed += l
                        pos_ref += l
            # attempt to get mismatches from MD tag
            md = safe_get_tag(read, 'MD')
            nm = safe_get_tag(read, 'NM')
            if md is not None:
                md_mismatches, md_del = parse_md_mismatches(md)
                mismatches = md_mismatches
                # md_del should equal del_bases
                # distribute mismatches across covered region proportionally: simple uniform distribution
                # but better: map MD to positions is complex; as approximation, distribute mismatches evenly across covered segments
                # We'll compute per-alignment distribution later: allocate mismatches to bins proportionally to covered lengths
                # For now store per-alignment counts, and distribute after
                # We will record covered regions to distribute mismatches below
                # collect covered segments for this read for later distribution
            else:
                # fallback: if NM present, compute mismatches = NM - indels
                if nm is not None:
                    try:
                        nmv = int(nm)
                        mismatches = max(0, nmv - (ins_bases + del_bases))
                    except:
                        mismatches = 0
                else:
                    mismatches = 0
            # For inaccurate distribution we will approx mismatches uniformly across covered blocks encountered
            # To do that, iterate again to collect covered segments and distribute mismatches
            # Collect covered segments of reference ("consuming ref", excluding D which we already counted)
            pos_ref = rstart0
            covered_blocks = []
            for op, l in cig:
                if op in (0,7,8):  # M or eqx
                    covered_blocks.append((pos_ref, pos_ref + l))
                    pos_ref += l
                elif op == 2:
                    pos_ref += l
                elif op == 1:
                    pass
                else:
                    if op in (3,):
                        pos_ref += l
            total_cov_len = sum(end - st for st, end in covered_blocks)
            if total_cov_len > 0 and mismatches > 0:
                for st, end in covered_blocks:
                    add_count_to_bins(chrom_bins[rname]['mismatch'], rlen[rname], bin_size, st, end, mismatches * ((end - st) / total_cov_len))
            # done with distribution for MD-case

        # note: we already added per-bin contributions for matches/mismatches/dels/ins during op parsing above
        # But in MD-case we added mismatches via add_count_to_bins; ensure matches accounting:
        # compute match = covered - mismatches - del
        # We'll compute matches per bin at end as covered - mismatch - del

        # For align-level summary compute covered/ref aln len
        aligned_ref_len = 0
        # aligned_ref_len should be sum of ops consuming reference: 0(M),2(D),7,8
        for op, l in cig:
            if op in (0,2,7,8):
                aligned_ref_len += l
        # For NM fallback we used nm
        # use NM tag as total edits if available
        nm_tag = safe_get_tag(read, 'NM')
        # sum ins/del already computed; mismatches may have been set from eqx or MD or NM fallback

        # Estimate matches:
        matches_est = max(0, aligned_ref_len - mismatches - del_bases)
        # add matches to bins (if eqx we already added matches; if not, add difference)
        if not has_eqx:
            # we previously added covered and mismatches/del in MD distribution,
            # but match bins not set - set matches = covered - mismatch - del via per-bin computation at the end
            pass

        # accumulate totals
        total_aligned += aligned_ref_len
        total_mism += mismatches
        total_ins += ins_bases
        total_del += del_bases

        # record alignment row
        pct_id = None
        denom = (mismatches + matches_est)
        if denom > 0:
            pct_id = matches_est / denom
        align_rows.append({
            'qname': qname, 'rname': rname,
            'rstart1': rstart0 + 1, 'rend1': rend0,
            'aligned_ref_len': aligned_ref_len,
            'mismatches': mismatches, 'ins_bases': ins_bases, 'del_bases': del_bases,
            'pct_identity': pct_id, 'NM_tag': nm_tag
        })

    # after processing reads, compute per-bin matches = covered - mismatch - del
    # and compute stats per bin
    windows = []
    genome_cov = 0.0
    genome_mism = 0.0
    genome_ins = 0.0
    genome_del = 0.0
    for chrom in refs:
        nb = n_bins[chrom]
        cov = chrom_bins[chrom]['covered']
        mism = chrom_bins[chrom]['mismatch']
        ins = chrom_bins[chrom]['ins']
        dele = chrom_bins[chrom]['del']
        # matches as covered - mismatch - del
        match = cov - mism - dele
        match[match < 0] = 0.0
        chrom_bins[chrom]['match'] = match
        for b in range(nb):
            start = b * bin_size + 1
            end = min(rlen[chrom], (b + 1) * bin_size)
            covered_b = float(cov[b])
            mism_b = float(mism[b])
            ins_b = float(ins[b])
            del_b = float(dele[b])
            divergence = mism_b / covered_b if covered_b > 0 else None
            indel_rate = (ins_b + del_b) / covered_b if covered_b > 0 else None
            windows.append({
                'chrom': chrom, 'start': start, 'end': end,
                'covered_bases': covered_b, 'mismatches': mism_b,
                'ins_bases': ins_b, 'del_bases': del_b,
                'divergence': divergence, 'indel_rate': indel_rate
            })
            genome_cov += covered_b
            genome_mism += mism_b
            genome_ins += ins_b
            genome_del += del_b

    # produce dataframes
    align_df = pd.DataFrame(align_rows)
    windows_df = pd.DataFrame(windows)
    genome_summary = pd.DataFrame([{
        'aligned_bases_total': genome_cov,
        'mismatches_total': genome_mism,
        'ins_total': genome_ins,
        'del_total': genome_del,
        'pairwise_pi': (genome_mism / genome_cov) if genome_cov > 0 else None,
        'indel_rate': (genome_ins + genome_del) / genome_cov if genome_cov > 0 else None
    }])

    # optional: aggregate by gene if GFF provided
    gene_df = None
    if gff_path:
        # read GFF lines, only keep features of type gene (or mRNA) columns: chrom, start, end, id
        feats = []
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                chrom_g, source, ftype, s, e, score, strand, phase, attr = parts[:9]
                if ftype.lower() not in ('gene','mrna','cds','exon', 'transcript'):
                    continue
                # try to extract ID or Name
                m = re.search(r'ID=([^;]+)', attr)
                gid = m.group(1) if m else (re.search(r'Name=([^;]+)', attr).group(1) if re.search(r'Name=([^;]+)', attr) else None)
                feats.append((chrom_g, int(s), int(e), gid if gid else f"{chrom_g}:{s}-{e}", ftype))
        if feats:
            # aggregate by gene by summing overlapping bins weighted by overlap
            rows = []
            # convert windows_df into per-chrom index for fast lookup
            win_group = windows_df.groupby('chrom')
            for chrom_g, s, e, gid, ftype in feats:
                if chrom_g not in win_group.groups:
                    continue
                dfw = win_group.get_group(chrom_g)
                # select bins overlapping feature
                overlap = dfw[(dfw['end'] >= s) & (dfw['start'] <= e)].copy()
                if overlap.empty:
                    continue
                # compute overlap length for weighting
                overlap['ovl'] = overlap.apply(lambda r: min(r['end'], e) - max(r['start'], s) + 1, axis=1)
                covered = (overlap['covered_bases'] * overlap['ovl'] / (overlap['end'] - overlap['start'] + 1)).sum()
                mism = (overlap['mismatches'] * overlap['ovl'] / (overlap['end'] - overlap['start'] + 1)).sum()
                insv = (overlap['ins_bases'] * overlap['ovl'] / (overlap['end'] - overlap['start'] + 1)).sum()
                delv = (overlap['del_bases'] * overlap['ovl'] / (overlap['end'] - overlap['start'] + 1)).sum()
                div = mism / covered if covered > 0 else None
                rows.append({'gene_id': gid, 'chrom': chrom_g, 'start': s, 'end': e,
                             'covered': covered, 'mismatches': mism, 'ins': insv, 'del': delv, 'divergence': div, 'type': ftype})
            gene_df = pd.DataFrame(rows)

    # write outputs
    align_df.to_csv(outprefix + ".alignments.csv", index=False)
    windows_df.to_csv(outprefix + ".windows.tsv", sep='\t', index=False)
    genome_summary.to_csv(outprefix + ".genome_summary.csv", index=False)
    if gene_df is not None:
        gene_df.to_csv(outprefix + ".genes.tsv", sep='\t', index=False)

    bam.close()

    print("[done] Outputs written with prefix:", outprefix)
    return

# ---------------------- CLI ----------------------
def main():
    parser = argparse.ArgumentParser(description="Compute haplotype divergence metrics from assembly-vs-assembly BAM")
    parser.add_argument("-b", "--bam", required=True, help="Sorted & indexed BAM from minimap2 (hap2->hap1)")
    parser.add_argument("-o", "--out", required=True, help="Output prefix")
    parser.add_argument("--bin", type=int, default=10000, help="Bin size for window summary (default 10000)")
    parser.add_argument("--min_mapq", type=int, default=20, help="Minimum alignment MAPQ to consider")
    parser.add_argument("--max_align", type=int, default=0, help="Max alignments to process (0=all)")
    parser.add_argument("--gff", default=None, help="Optional gene GFF/GTF to aggregate divergence by gene")
    args = parser.parse_args()

    max_align = None if args.max_align == 0 else args.max_align
    process_bam(args.bam, args.out, bin_size=args.bin, min_mapq=args.min_mapq, max_align_records=max_align, gff_path=args.gff)

if __name__ == "__main__":
    main()