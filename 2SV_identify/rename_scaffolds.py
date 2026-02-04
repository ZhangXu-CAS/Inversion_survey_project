#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import subprocess

# ========================
# Parameter settings
# ========================
if len(sys.argv) != 2:
    print("Usage: python rename_scaffolds.py <sample_id>")
    sys.exit(1)

SAMPLE_ID = sys.argv[1]

# Input files and directories
BEST_MATCHES_FILE = "best_matches.tsv"
HAP1_FOLDER = "hap1_scaffolds"
HAP2_FOLDER = "hap2_scaffolds"

# Output files
HAP1_OUT_FILE = f"{SAMPLE_ID}_hap1_renchr.fa"
HAP2_OUT_FILE = f"{SAMPLE_ID}_hap2_renchr.fa"

# ========================
# Step 1: Load and filter match table
# ========================
print("Reading and filtering best_matches.tsv...")
df = pd.read_csv(BEST_MATCHES_FILE, sep="\t", header=None, names=["hap1", "hap2", "distance"])
df = df[df["distance"] <= 0.1].copy()

# Remove hap2 duplicates to ensure one-to-one mapping
df = df.drop_duplicates(subset=["hap2"], keep="first")

# ========================
# Step 2: Extract hap1 scaffolds (in given order)
# ========================
print("Reading hap1 scaffolds...")
hap1_records = []
hap1_names = []

for hap1_file in df["hap1"]:
    path = os.path.join(HAP1_FOLDER, os.path.basename(hap1_file))
    for record in SeqIO.parse(path, "fasta"):
        hap1_records.append(record)
        hap1_names.append(record.id)

SeqIO.write(hap1_records, HAP1_OUT_FILE, "fasta")
print(f"Wrote {HAP1_OUT_FILE} with {len(hap1_records)} sequences.")

# ========================
# Step 3: Extract and rename hap2 scaffolds to match hap1 names
# ========================
print("Reading and renaming hap2 scaffolds...")
hap2_records = []

for hap1_file, hap2_file in zip(df["hap1"], df["hap2"]):
    new_name = os.path.basename(hap1_file).replace(".fa", "")
    hap2_path = os.path.join(HAP2_FOLDER, os.path.basename(hap2_file))
    
    for record in SeqIO.parse(hap2_path, "fasta"):
        record.id = new_name
        record.description = ""
        hap2_records.append(record)

if len(hap2_records) != len(hap1_records):
    print(f"Error: mismatched counts: hap1 = {len(hap1_records)}, hap2 = {len(hap2_records)}")
    sys.exit(1)

SeqIO.write(hap2_records, HAP2_OUT_FILE, "fasta")
print(f"Wrote {HAP2_OUT_FILE} with {len(hap2_records)} sequences.")

# ========================
# Step 4: Index output FASTAs
# ========================
print("Indexing fasta files...")
subprocess.run(f"samtools faidx {HAP1_OUT_FILE}", shell=True, check=True)
subprocess.run(f"samtools faidx {HAP2_OUT_FILE}", shell=True, check=True)

print("All done!")
