#!/bin/bash
# pipeline for SV detection using minimap + syri -- Xu Zhang 04/01/2025--add /tmp for parallel

# Exit immediately if any command exits with a non-zero status
set -e
trap 'echo "❗️Error: $(date) - $0 - line $LINENO" | tee -a $LOG_FILE' ERR

# Check for correct input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: bash pipe_minimap_syri.sh <sample_id>"
    exit 1
fi

SAMPLE_ID=$1
WD=`pwd`
mkdir -p $WD/tmp
cd $WD/$SAMPLE_ID
LOG_FILE="pipe_minimap_syri.log"

# Create necessary output directories
mkdir -p minimap_results syri_results hap1_scaffolds hap2_scaffolds hap1_mash_sketches hap2_mash_sketches syri_results_reminimap

echo "🚀 Starting structural variant pipeline for ${SAMPLE_ID}" | tee -a $LOG_FILE

# Step 1: Split hap1 and hap2 into individual scaffold FASTA files
echo "📚 Splitting FASTA files into scaffolds..." | tee -a $LOG_FILE

# Run Python script to split FASTA files into individual scaffold files
python3 - <<EOF
from Bio import SeqIO
import os

def split_fasta(input_fasta, output_dir):
    """Splits a FASTA file into individual scaffold files"""
    os.makedirs(output_dir, exist_ok=True)
    for record in SeqIO.parse(input_fasta, "fasta"):
        out_file = os.path.join(output_dir, f"{record.id}.fa")
        SeqIO.write(record, out_file, "fasta")

split_fasta("${SAMPLE_ID}_hap1_chr.fa", "hap1_scaffolds")
split_fasta("${SAMPLE_ID}_hap2_chr.fa", "hap2_scaffolds")
EOF

echo "✅ FASTA splitting completed." | tee -a $LOG_FILE

# Step 2: Generate Mash sketches in parallel
echo "⚡️ Generating Mash sketches..." | tee -a $LOG_FILE

# Use GNU parallel to accelerate mash sketch creation
find hap1_scaffolds -name "*.fa" | parallel --tmpdir $WD/tmp -j 20 "mash sketch -o hap1_mash_sketches/{/.} {}"
find hap2_scaffolds -name "*.fa" | parallel --tmpdir $WD/tmp -j 20 "mash sketch -o hap2_mash_sketches/{/.} {}"

# Step 3: Compute Mash distances
echo "🔎 Computing Mash distances..." | tee -a $LOG_FILE
MASH_OUT="mash_distances.txt"
rm -f $MASH_OUT

# Use parallel to compute Mash distances between hap1 and hap2 sketches
parallel  --tmpdir $WD/tmp -j 20 "mash dist {} hap2_mash_sketches/*.msh >> $MASH_OUT" ::: hap1_mash_sketches/*.msh

# Step 4: Identify homologous scaffold pairs using Hungarian algorithm
echo "🧩 Identifying homologous scaffold pairs..." | tee -a $LOG_FILE
python3 - <<EOF
from collections import defaultdict
def find_min_matches(input_file):
    data = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('[file name]') or line.startswith('[file content'):
                continue 
            parts = line.split('\t')
            if len(parts) < 3:
                continue 
            hap1 = parts[0]
            hap2 = parts[1]
            try:
                distance = float(parts[2])
            except ValueError:
                continue 
            data.append((hap1, hap2, distance))
    groups = defaultdict(list)
    for entry in data:
        groups[entry[0]].append(entry)
    results = []
    for hap1, entries in groups.items():
        min_entry = min(entries, key=lambda x: x[2])
        results.append((hap1, min_entry[1], min_entry[2]))
    # sort hap1 scaffold
    results.sort(key=lambda x: x[0])
    return results
def main():
    input_file = "mash_distances.txt"
    output_file = "best_matches.tsv"
    matches = find_min_matches(input_file)
    with open(output_file, 'w') as f:
        for hap1, hap2, dist in matches:
            f.write(f"{hap1}\t{hap2}\t{dist}\n")
if __name__ == "__main__":
    main()
EOF

echo "✅ Homologous scaffold pairs identified." | tee -a $LOG_FILE

# Step 5: Rename hap2 scaffolds to match hap1

echo "🔄 Renaming and sorting hap2 scaffolds..." | tee -a $LOG_FILE
python3 $WD/rename_scaffolds.py ${SAMPLE_ID}

# Step 6: Run first-round minimap + syri
echo "⚔️ Running first minimap2 + syri..." | tee -a $LOG_FILE
minimap2 -t 40 -ax asm5 --eqx ${SAMPLE_ID}_hap1_renchr.fa ${SAMPLE_ID}_hap2_renchr.fa > minimap_results/1st_chr_alignment.sam 

# Activate syri environment
source  /home/xuzhang/software/miniconda3/bin/activate syri
ulimit -v 200000000
syri -c minimap_results/1st_chr_alignment.sam \
     -r ${SAMPLE_ID}_hap1_renchr.fa \
     -q ${SAMPLE_ID}_hap2_renchr.fa \
     --dir syri_results \
     -k -F S 2> syri_results/syri_error_out.txt

# Extract scaffolds with high inversion fractions for reversal
cat syri_results/syri_error_out.txt | grep WARN | grep "Reference chromosome" | cut -f 23 -d " " | sed s/\(//g | sed s/\)\.//g > chr_to_rev.txt

# Step 7: Reverse hap1 scaffolds and redo minimap + syri
echo "🔁 Reversing and re-aligning hap1 scaffolds..." | tee -a $LOG_FILE
python3 $WD/reverse_scaffolds.py -i ${SAMPLE_ID}_hap1_renchr.fa -r chr_to_rev.txt

minimap2 -t 40 -ax asm5 --eqx ${SAMPLE_ID}_hap1_revchr.fa ${SAMPLE_ID}_hap2_renchr.fa > minimap_results/2nd_chr_alignment.sam 

syri -c minimap_results/2nd_chr_alignment.sam \
     -r ${SAMPLE_ID}_hap1_revchr.fa \
     -q ${SAMPLE_ID}_hap2_renchr.fa \
     --dir syri_results_reminimap \
     -k -F S 2> syri_results_reminimap/syri_error_out.txt

# Step 8: Plot results for visualization
echo "📊 Generating plot for genome alignment..." | tee -a $LOG_FILE
echo -e "${SAMPLE_ID}_hap1_revchr.fa\t${SAMPLE_ID}_hap1\n${SAMPLE_ID}_hap2_renchr.fa\t${SAMPLE_ID}_hap2" > genomes.txt
plotsr --sr syri_results_reminimap/syri.out --genomes genomes.txt -o ${SAMPLE_ID}_minimap_syri_plot.png

# Step 9: Run nucmer + syri for fine resolution
# mkdir -p chr_numcmer syri_chr_numcmer
#echo "🧬 Running nucmer + syri for final alignment..." | tee -a $LOG_FILE
#nucmber_out="chr_numcmer/rechr"
#
## Run nucmer and filter alignments
#nucmer --maxmatch -c 100 -b 500 -l 50 ${SAMPLE_ID}_hap1_revchr.fa ${SAMPLE_ID}_hap2_renchr.fa -p $nucmber_out
#delta-filter -m -i 90 -l 100 $nucmber_out.delta > $nucmber_out.filtered.delta
#show-coords -THrd $nucmber_out.filtered.delta > $nucmber_out.filtered.coords
#
## Run syri for the final refined analysis
#syri -c $nucmber_out.filtered.coords -d $nucmber_out.filtered.delta \
#     -r ${SAMPLE_ID}_hap1_revchr.fa \
#     -q ${SAMPLE_ID}_hap2_renchr.fa \
#     --dir syri_chr_numcmer -k --nc 12 2> syri_chr_numcmer/syri_error_out.txt
source  /home/xuzhang/miniconda3/bin/deactivate 
echo "✅ All steps completed successfully! Results are available in the respective directories." | tee -a $LOG_FILE
