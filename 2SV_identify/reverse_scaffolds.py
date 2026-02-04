#!/usr/bin/env python3
"""
reverse_scaffolds.py

This script reverses and complements selected scaffolds in a hap1 FASTA file.

Usage:
    python reverse_scaffolds.py -i <hap1_chr.fa> -r chr_to_rev.txt
"""

import argparse
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_reverse_list(reverse_file):
    """Load list of scaffolds to reverse from the file."""
    with open(reverse_file, 'r') as f:
        return {line.strip() for line in f}

def process_fasta(input_file, output_file, reverse_set):
    """Reverse complementary of selected scaffolds and write to new FASTA."""
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        scaffold_id = record.id.split()[0]
        if scaffold_id in reverse_set:
            reversed_seq = record.seq.reverse_complement()
            new_record = SeqRecord(
                reversed_seq,
                id=record.id,
                description=record.description + " [revcomp]"
            )
        else:
            new_record = record
        records.append(new_record)

    with open(output_file, 'w') as f:
        SeqIO.write(records, f, "fasta")
        
def index_fasta(output_file):
    """Create an index file for the output FASTA using samtools."""
    subprocess.run(['samtools', 'faidx', output_file], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reverse complement specific scaffolds in FASTA.")
    parser.add_argument('-i', '--input', required=True,
                        help="Input hap1 genome FASTA file (format: $id_hap1_renchr.fa)")
    parser.add_argument('-r', '--reverse', default='chr_to_rev.txt',
                        help="List of scaffolds to reverse (default: chr_to_rev.txt)")
    
    args = parser.parse_args()

    base_name = os.path.basename(args.input)
    if '_hap1_renchr.fa' not in base_name:
        raise ValueError("Input file must match the format: $id_hap1_chr.fa")
    
    output_file = base_name.replace('_hap1_renchr.fa', '_hap1_revchr.fa')

    print("Loading reverse complement list...")
    reverse_set = load_reverse_list(args.reverse)

    print(f"Processing {args.input} -> {output_file}...")
    process_fasta(args.input, output_file, reverse_set)

    print(f"✅ Done! Output: {output_file}")