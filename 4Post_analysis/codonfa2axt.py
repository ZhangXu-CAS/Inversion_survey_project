#!/usr/bin/env python3
"""
codonfa2axt.py
Convert a two-sequence codon-aligned FASTA to AXT format (single block).
Usage:
    python3 codonfa2axt.py in.codon.fa out.axt [query_name] [subject_name]

Requirements:
    biopython
Notes:
    - Input FASTA must contain exactly two sequences (codon-aligned; gaps as '-').
    - The script will remove gaps for coordinate calculation but keep alignment in AXT.
    - AXT format: <index> <qname> <qstart> <qend> <sname> <sstart> <send> <strand> <score>
                <qseq>
                <sseq>
    - We set index=1 and strand='+' and score=0 by default.
"""
import sys
from Bio import SeqIO
import re

if len(sys.argv) < 3:
    print("Usage: codonfa2axt.py in.codon.fa out.axt [qname] [sname]", file=sys.stderr)
    sys.exit(1)

inf = sys.argv[1]
outf = sys.argv[2]
qname_arg = sys.argv[3] if len(sys.argv) > 3 else None
sname_arg = sys.argv[4] if len(sys.argv) > 4 else None

records = list(SeqIO.parse(inf, "fasta"))
if len(records) != 2:
    print("Error: input fasta must contain exactly 2 sequences.", file=sys.stderr)
    sys.exit(1)

qrec = records[0]
srec = records[1]
qname = qname_arg if qname_arg else qrec.id.split()[0]
sname = sname_arg if sname_arg else srec.id.split()[0]

qseq = str(qrec.seq).replace('\n','').strip()
sseq = str(srec.seq).replace('\n','').strip()

# sanity: lengths equal and divisible by 3 (not strictly required for AXT but good)
if len(qseq) != len(sseq):
    print("Warning: aligned sequences have different lengths.", file=sys.stderr)

# compute ungapped coordinates (1-based inclusive)
def ungapped_coords(aln_seq):
    # remove gaps '-' and count nucleotides -> positions in original CDS
    ungapped = aln_seq.replace('-', '')
    # start is 1, end is len(ungapped)
    return 1, max(1, len(ungapped))

qstart, qend = ungapped_coords(qseq)
sstart, send = ungapped_coords(sseq)

# Remove trailing newlines from sequences and ensure no lowercase
qseq_line = qseq
sseq_line = sseq

# Write AXT: a single block, index=1, strand='+', score=0
with open(outf, 'w') as fo:
    fo.write(f"1 {qname} {qstart} {qend} {sname} {sstart} {send} + 0\n")
    fo.write(qseq_line + "\n")
    fo.write(sseq_line + "\n\n")
print("Wrote", outf)
