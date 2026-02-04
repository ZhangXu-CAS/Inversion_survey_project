#!/usr/bin/env python3
"""
pav_from_blast.py

Usage:
  python3 pav_from_blast.py \
    --hap1-prot hap1.prot.fa --hap2-prot hap2.prot.fa \
    --b1 hap1_vs_hap2.blast --b2 hap2_vs_hap1.blast \
    --outdir pav_blast_out 
   

Outputs in outdir:
  - pav_summary.tsv (hap1_total, hap2_total, rbh_pairs, hap1_only, hap2_only, 1:many, many:1, many:many)
  - per_gene_presence_hap1.tsv (for each hap1 gene: copy_in_hap2, present_flag, top_hits...)
  - per_gene_presence_hap2.tsv
  - rbh_pairs.tsv (hap1_id \t hap2_id)
  - pair_classification.tsv (hap1_id, hap2_id, class)
"""
import argparse, os, sys, csv
from collections import defaultdict

def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--hap1-prot', required=True)
    p.add_argument('--hap2-prot', required=True)
    p.add_argument('--b1', required=True, help='hap1_vs_hap2.blast (q=hap1)')
    p.add_argument('--b2', required=True, help='hap2_vs_hap1.blast (q=hap2)')
    p.add_argument('--outdir', default='pav_blast_out')
    p.add_argument('--pid', type=float, default=20.0)
    p.add_argument('--cov_q', type=float, default=0.5)
    p.add_argument('--cov_s', type=float, default=0.5)
    p.add_argument('--e', type=float, default=1e-3)
    return p.parse_args()

def read_fasta_ids_lengths(fa):
    ids = []
    lengths = {}
    with open(fa) as fh:
        name=None; seq=[]
        for line in fh:
            if line.startswith('>'):
                if name:
                    s=''.join(seq)
                    ids.append(name)
                    lengths[name]=len(s)
                name=line[1:].strip().split()[0]
                seq=[]
            else:
                seq.append(line.strip())
        if name:
            s=''.join(seq)
            ids.append(name)
            lengths[name]=len(s)
    return ids, lengths

def parse_blast_tab(fn, pid, cov_q, cov_s, emax):
    """
    Expect blast outfmt 6 with fields:
    qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
    Return dict: hits[q].append((s,pident,aln_len,qlen,slen,evalue,bitscore))
    """
    hits=defaultdict(list)
    with open(fn) as fh:
        for line in fh:
            parts=line.rstrip().split('\t')
            if len(parts) < 12: 
                continue
            q,s = parts[0], parts[1]
            pident = float(parts[2])
            aln_len = int(parts[3])
            qlen = int(parts[4])
            slen = int(parts[5])
            e = float(parts[10])
            bitscore = float(parts[11])
            # coverage checks
            covq = aln_len / qlen if qlen>0 else 0
            covs = aln_len / slen if slen>0 else 0
            if pident >= pid and covq >= cov_q and covs >= cov_s and e <= emax:
                hits[q].append((s,pident,aln_len,qlen,slen,e,bitscore))
    return hits

def get_top_hit_map(hits):
    # return q->top subject (highest bitscore)
    top={}
    for q,hs in hits.items():
        if not hs: continue
        top_hit = max(hs, key=lambda x: (x[-1], x[1])) # by bitscore then pid
        top[q]=top_hit[0]
    return top

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    hap1_ids, hap1_len = read_fasta_ids_lengths(args.hap1_prot)
    hap2_ids, hap2_len = read_fasta_ids_lengths(args.hap2_prot)

    hap1_set = set(hap1_ids); hap2_set=set(hap2_ids)
    hap1_total = len(hap1_ids); hap2_total=len(hap2_ids)

    hits_h1_to_h2 = parse_blast_tab(args.b1, args.pid, args.cov_q, args.cov_s, args.e)
    hits_h2_to_h1 = parse_blast_tab(args.b2, args.pid, args.cov_q, args.cov_s, args.e)

    # presence: any passing hit means present
    hap1_present_in_hap2 = {}
    for gid in hap1_ids:
        hap1_present_in_hap2[gid] = sorted([h for h in hits_h1_to_h2.get(gid,[])], key=lambda x: (-x[-1], -x[1]))
    hap2_present_in_hap1 = {}
    for gid in hap2_ids:
        hap2_present_in_hap1[gid] = sorted([h for h in hits_h2_to_h1.get(gid,[])], key=lambda x: (-x[-1], -x[1]))

    # top-hit maps for RBH detection
    top_h1 = get_top_hit_map(hits_h1_to_h2)
    top_h2 = get_top_hit_map(hits_h2_to_h1)

    # find RBH: hap1->s and s->hap1
    rbh_pairs=[]
    for q, s in top_h1.items():
        if s in top_h2 and top_h2[s] == q:
            rbh_pairs.append((q,s))

    # classify each RBH pair as 1:1 or part of many-to-many later
    # For copy counts: count number of subject hits passing thresholds
    per_hap1=[]
    hap1_only_count=0
    for gid in hap1_ids:
        hits = hap1_present_in_hap2.get(gid,[])
        copy_num = len({h[0] for h in hits})
        present_flag = 1 if copy_num>0 else 0
        if present_flag==0:
            hap1_only_count += 1
        top_hits_str = ",".join([f"{h[0]}|pid={h[1]:.1f}|aln={h[2]}" for h in hits])
        per_hap1.append((gid, copy_num, present_flag, top_hits_str))

    per_hap2=[]
    hap2_only_count=0
    for gid in hap2_ids:
        hits = hap2_present_in_hap1.get(gid,[])
        copy_num = len({h[0] for h in hits})
        present_flag = 1 if copy_num>0 else 0
        if present_flag==0:
            hap2_only_count += 1
        top_hits_str = ",".join([f"{h[0]}|pid={h[1]:.1f}|aln={h[2]}" for h in hits])
        per_hap2.append((gid, copy_num, present_flag, top_hits_str))

    # classify pair relationship categories (1:1, 1:many, many:1, many:many)
    # Generate maps:
    h1_to_subjects = {gid: set([h[0] for h in hits_h1_to_h2.get(gid,[])]) for gid in hap1_ids}
    h2_to_subjects = {gid: set([h[0] for h in hits_h2_to_h1.get(gid,[])]) for gid in hap2_ids}
    # invert hap1->subject to subject->set(h1)
    subj_to_h1 = defaultdict(set)
    for h1,gset in h1_to_subjects.items():
        for s in gset:
            subj_to_h1[s].add(h1)
    # invert hap2->subject (subject is hap1 id)
    subj_to_h2 = defaultdict(set)
    for h2,gset in h2_to_subjects.items():
        for s in gset:
            subj_to_h2[s].add(h2)

    # classify RBH-like groups by checking sizes
    one2one=0; one2many=0; many2one=0; many2many=0
    pair_rows=[]
    # iterate through RBH pairs first
    seen_pairs=set()
    for h1,h2 in rbh_pairs:
        # number of hap2 hits for h1:
        n_h2 = len(h1_to_subjects.get(h1,[]))
        # number of hap1 hits pointing to h2:
        n_h1 = len(subj_to_h1.get(h2,[]))
        if n_h2==1 and n_h1==1:
            one2one += 1
            cls='1:1'
        elif n_h2>1 and n_h1==1:
            one2many += 1
            cls='1:many'
        elif n_h2==1 and n_h1>1:
            many2one += 1
            cls='many:1'
        else:
            many2many += 1
            cls='many:many'
        pair_rows.append((h1,h2,cls,n_h2,n_h1))
        seen_pairs.add((h1,h2))

    # Also consider non-RBH relationships if you want; here we focus on RBH classification

    # write outputs
    outdir=args.outdir
    with open(os.path.join(outdir,'per_gene_presence_hap1.tsv'),'w') as fo:
        fo.write("hap1_id\tcopy_in_hap2\tpresent_flag\ttop_hits\n")
        for row in per_hap1:
            fo.write("\t".join(map(str,row)) + "\n")
    with open(os.path.join(outdir,'per_gene_presence_hap2.tsv'),'w') as fo:
        fo.write("hap2_id\tcopy_in_hap1\tpresent_flag\ttop_hits\n")
        for row in per_hap2:
            fo.write("\t".join(map(str,row)) + "\n")
    with open(os.path.join(outdir,'rbh_pairs.tsv'),'w') as fo:
        for a,b in rbh_pairs:
            fo.write(f"{a}\t{b}\n")
    with open(os.path.join(outdir,'pair_classification.tsv'),'w') as fo:
        fo.write("hap1_id\thap2_id\tclass\th1->h2_count\th2->h1_count\n")
        for r in pair_rows:
            fo.write("\t".join(map(str,r)) + "\n")

    with open(os.path.join(outdir,'pav_summary.tsv'),'w') as fo:
        fo.write("metric\tcount\n")
        fo.write(f"hap1_total\t{hap1_total}\n")
        fo.write(f"hap2_total\t{hap2_total}\n")
        fo.write(f"rbh_pairs\t{len(rbh_pairs)}\n")
        fo.write(f"hap1_only\t{hap1_only_count}\n")
        fo.write(f"hap2_only\t{hap2_only_count}\n")
        fo.write(f"1:1\t{one2one}\n")
        fo.write(f"1:many\t{one2many}\n")
        fo.write(f"many:1\t{many2one}\n")
        fo.write(f"many:many\t{many2many}\n")

    print("Wrote outputs to", outdir)
    print("Summary: hap1_total",hap1_total,"hap2_total",hap2_total,"rbh",len(rbh_pairs),
          "hap1_only",hap1_only_count,"hap2_only",hap2_only_count)

if __name__=='__main__':
    main()