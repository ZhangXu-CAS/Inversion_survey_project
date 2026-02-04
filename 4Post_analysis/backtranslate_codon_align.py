#!/usr/bin/env python3
"""
backtranslate_codon_align.py

用法:
python3 backtranslate_codon_align.py protein_aln.fa hap1_cds.fa hap2_cds.fa id1 id2 out_codon_aln.fa

要求:
- protein_aln.fa: 两条蛋白序列的多序列比对文件 (fasta)，包含 id1 和 id2
- hap?_cds.fa: 原始 CDS nucleotide fasta，header 与蛋白 ID 一致（或者能通过前缀匹配）
输出:
- codon alignment fasta (核酸，三联码对齐)
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq

if len(sys.argv) != 7:
    print("Usage: backtranslate_codon_align.py prot_aln.fa hap1_cds.fa hap2_cds.fa id1 id2 out.fa")
    sys.exit(1)

prot_aln_fa, hap1_cds_fa, hap2_cds_fa, id1, id2, outfa = sys.argv[1:]

# load CDS dicts
cds1 = {rec.id: str(rec.seq) for rec in SeqIO.parse(hap1_cds_fa, "fasta")}
cds2 = {rec.id: str(rec.seq) for rec in SeqIO.parse(hap2_cds_fa, "fasta")}

# load protein alignment
prot_aln = {rec.id: str(rec.seq) for rec in SeqIO.parse(prot_aln_fa, "fasta")}

if id1 not in prot_aln or id2 not in prot_aln:
    # try header prefix matches
    keys = list(prot_aln.keys())
    # fallback: look for entries starting with id
    found1 = next((k for k in keys if k.split()[0]==id1 or k.startswith(id1)), None)
    found2 = next((k for k in keys if k.split()[0]==id2 or k.startswith(id2)), None)
    if found1 and found2:
        id1 = found1
        id2 = found2
    else:
        print("Error: ids not found in protein alignment:", id1, id2, file=sys.stderr)
        sys.exit(1)

p1 = prot_aln[id1]
p2 = prot_aln[id2]

# find nucleotide sequences for original proteins
# assume protein id maps directly to cds id; if not, try prefix match
def find_cds(cds_dict, qid):
    if qid in cds_dict:
        return cds_dict[qid]
    # try prefix match
    for k in cds_dict:
        if k.split()[0] == qid or k.startswith(qid):
            return cds_dict[k]
    return None

n1 = find_cds(cds1, id1)
n2 = find_cds(cds2, id2)

if n1 is None or n2 is None:
    print("Error: cannot find CDS for", id1, id2, file=sys.stderr)
    sys.exit(1)

# translate CDS to protein (without stop)
prot_from_cds1 = str(Seq(n1).translate(to_stop=False))
prot_from_cds2 = str(Seq(n2).translate(to_stop=False))

# Now, build codon alignment by iterating aligned protein sequences and mapping back to codons
def backtranslate(aligned_prot, cds_nt):
    codon_list = []
    cds_index = 0
    # remove potential trailing stop codon in cds translation? we assume cds length divisible by 3
    for aa in aligned_prot:
        if aa == '-':
            codon_list.append('---')
        else:
            codon = cds_nt[cds_index:cds_index+3]
            if len(codon) < 3:
                # padding
                codon = codon + '-'*(3-len(codon))
            codon_list.append(codon)
            cds_index += 3
    return ''.join(codon_list)

codon_seq1 = backtranslate(p1, n1)
codon_seq2 = backtranslate(p2, n2)

# sanity check: lengths equal and divisible by 3?
if len(codon_seq1) != len(codon_seq2):
    print("Warning: codon seqs different length", file=sys.stderr)

# write fasta
with open(outfa, 'w') as fh:
    fh.write(f">{id1}\n")
    # wrap 60
    for i in range(0, len(codon_seq1), 60):
        fh.write(codon_seq1[i:i+60] + "\n")
    fh.write(f">{id2}\n")
    for i in range(0, len(codon_seq2), 60):
        fh.write(codon_seq2[i:i+60] + "\n")