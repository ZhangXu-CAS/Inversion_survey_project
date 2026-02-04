#!/bin/bash
# Script for HiFi + Hi-C assembly using hifiasm

set -euo pipefail  # Enable strict error handling

# Set working directory
wd="/Volumes/70T/AsteroidScratch/xuzhang/inver_proj/0data/ncbi"
cd "$wd" || { echo "Working directory $wd does not exist."; exit 1; }

# Ensure the 'ids' file exists
ids_file="ids"
if [[ ! -f "$ids_file" ]]; then
    echo "IDs file '$ids_file' not found in $wd."
    exit 1
fi

# Check if required tools are installed
required_tools=("hifiasm" "gfatools" "samtools" "seqkit")
for tool in "${required_tools[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool is not installed or not in PATH."
        exit 1
    fi
done

# Process each ID
while IFS= read -r id || [[ -n "$id" ]]; do
    echo "Processing ID: $id at $(date)"
    
    project_dir="$wd/$id"
    
    # Check if project directory exists
    if [[ ! -d "$project_dir" ]]; then
        echo "Directory $project_dir does not exist. Skipping."
        continue
    fi
    
    cd "$project_dir" || { echo "Failed to navigate to $project_dir."; continue; }
    
    # Create output directories if they don't exist
    mkdir -p hifiasm_hic_out asm_hifiasm_hic
    
    # Run hifiasm if the output does not already exist
    hifiasm_output="hifiasm_hic_out/${id}.hifiasm"
    hifiasm_log="${id}.hifiasm_hic.log"
    THREADS=40
	
    if [[ ! -f "$hifiasm_output.ctg.gfa" ]]; then
        echo "Running hifiasm for $id..."
        hifiasm -o "$hifiasm_output" -t "$THREADS"  \
            --h1 hic/*_1.fastq.gz \
            --h2 hic/*_2.fastq.gz \
            hifi/*.fastq.gz 2> "$hifiasm_log"
    else
        echo "hifiasm output for $id already exists. Skipping hifiasm step."
    fi
    
    # Convert GFA to FASTA for hap1
    hap1_gfa="hifiasm_hic_out/*.hap1.p_ctg.gfa"
    hap1_fa="asm_hifiasm_hic/${id}.hifiasm.hic.hap1.fa"
    if [[ ! -f "$hap1_fa" ]]; then
        echo "Converting GFA to FASTA for hap1 of $id..."
        gfatools gfa2fa $hap1_gfa > "$hap1_fa"
        samtools faidx "$hap1_fa"
    else
        echo "FASTA for hap1 of $id already exists. Skipping conversion."
    fi
    
    # Convert GFA to FASTA for hap2
    hap2_gfa="hifiasm_hic_out/*.hap2.p_ctg.gfa"
    hap2_fa="asm_hifiasm_hic/${id}.hifiasm.hic.hap2.fa"
    if [[ ! -f "$hap2_fa" ]]; then
        echo "Converting GFA to FASTA for hap2 of $id..."
        gfatools gfa2fa $hap2_gfa > "$hap2_fa"
        samtools faidx "$hap2_fa"
    else
        echo "FASTA for hap2 of $id already exists. Skipping conversion."
    fi
    
    # Generate statistics using seqkit
    stats_file="asm_hifiasm_hic/2hap.stat"
    if [[ ! -f "$stats_file" ]]; then
        echo "Generating statistics for $id..."
        seqkit stat -Ta "$hap1_fa" "$hap2_fa" > "$stats_file"
    else
        echo "Statistics for $id already exist. Skipping seqkit stat."
    fi
    
    # Move GFA files to the project directory if not already moved
    if ls hifiasm_hic_out/*.hap*.p_ctg.gfa 1> /dev/null 2>&1; then
        echo "Moving GFA files for $id..."
        mv hifiasm_hic_out/*.hap*.p_ctg.gfa ./
    else
        echo "No GFA files to move for $id."
    fi
    
    # Optional: Compress the hifiasm_hic_out directory
    # Uncomment the following lines if you want to archive the directory
    # echo "Compressing hifiasm_hic_out for $id..."
    # tar -cf - hifiasm_hic_out | pigz -p 50 > hifiasm_hic_out.tar.gz
    
    # Clean up the hifiasm_hic_out directory
    echo "Removing hifiasm_hic_out directory for $id..."
    rm -rf hifiasm_hic_out
    
    echo "$id hifiasm+hic done at $(date)"
    
    # Return to the working directory
    cd "$wd" || { echo "Failed to navigate back to $wd."; exit 1; }
    
done < "$ids_file"

echo "All projects processed successfully."
