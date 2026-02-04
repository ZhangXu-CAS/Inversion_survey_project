#!/usr/bin/env python3
"""
merge_codonfas_to_axt.py

Convert all two-sequence codon-aligned FASTA files in a directory into a single merged AXT file
with incrementing block indices.

Usage:
    python3 merge_codonfas_to_axt.py codon_dir merged.axt

Assumptions:
 - Each input file contains exactly two sequences (aligned codons with '-' for gaps).
 - Filename contains the two IDs separated by '_' (e.g. id1_id2.codon.aln.fa). If not,
   the script will use fasta headers (split by whitespace) as qname/sname.
 - The script writes AXT blocks with simplistic coordinates: qstart=1, qend=ungapped_len, similarly for subject.
"""
import sys, os
from pathlib import Path
from Bio import SeqIO

def ungapped_len(seq):
    return len(seq.replace('-', '').replace('\n',''))

def main():
    if len(sys.argv) != 3:
        print("Usage: merge_codonfas_to_axt.py <codon_dir> <out_axt>", file=sys.stderr)
        sys.exit(1)
    codon_dir = Path(sys.argv[1])
    out_axt = Path(sys.argv[2])
    files = sorted(codon_dir.glob("*codon*.fa*"))  # match many patterns

    idx = 1
    with out_axt.open('w') as fo:
        for f in files:
            recs = list(SeqIO.parse(str(f), "fasta"))
            if len(recs) != 2:
                print(f"Skipping {f}: expected 2 sequences, found {len(recs)}", file=sys.stderr)
                continue
            # try to extract qname and sname from filename if possible (id1__id2...)
            fname = f.name
            if "_" in fname:
                qname,sname = fname.split("_",1)[0], fname.split("_",1)[1].split('.')[0]
            else:
                qname = recs[0].id.split()[0]
                sname = recs[1].id.split()[0]

            qseq = str(recs[0].seq).replace('\r','')
            sseq = str(recs[1].seq).replace('\r','')

            q_ungap = ungapped_len(qseq)
            s_ungap = ungapped_len(sseq)
            # if ungapped lengths are zero skip
            if q_ungap == 0 or s_ungap == 0:
                print(f"Skipping {f}: ungapped length zero for one sequence", file=sys.stderr)
                continue

            # AXT header: idx qname qstart qend sname sstart send strand score
            # We'll use qstart=1, qend=q_ungap, sstart=1, send=s_ungap, strand='+', score=0
            fo.write(f"{idx} {qname} {1} {q_ungap} {sname} {1} {s_ungap} + 0\n")
            fo.write(qseq + "\n")
            fo.write(sseq + "\n\n")
            idx += 1

    print(f"Wrote {out_axt} with {idx-1} blocks")

if __name__ == "__main__":
    main()