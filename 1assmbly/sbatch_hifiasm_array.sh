#!/bin/sh
#SBATCH --job-name=1227hifi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=40:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-59  


# Read the ID for this task from the file 'ids'
ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ids)

# Navigate to the working directory
WORK_DIR="/lustre07/scratch/zhangxu/inversion_proj/1212_updates/$ID"
cd "$WORK_DIR" || { echo "Error: Directory $WORK_DIR not found!" >> "$ID.error.log"; exit 1; }

# Create output directories if they do not exist
mkdir -p hifiasm_hic_out asm_hifiasm_hic

echo "Starting run for $ID at: $(date)" >> "$ID.hifi.log"

# Run hifiasm with specified parameters
echo "Running hifiasm for $ID..." >> "$ID.hifi.log"
hifiasm -o hifiasm_hic_out/"$ID".hifiasm -t 40 --h1 hic_data/*_1.fastq.gz --h2 hic_data/*_2.fastq.gz hifi_data/*.fastq.gz 2>>"$ID.hifi.log"
echo "$ID hifiasm completed at $(date)" >> "$ID.hifi.log"

# Generate statistics for the assemblies
echo "Generating statistics for $ID..." >> "$ID.hifi.log"
gfatools gfa2fa hifiasm_hic_out/*.hap1.p_ctg.gfa > asm_hifiasm_hic/$ID.hifiasm.hic.hap1.fa 
samtools faidx asm_hifiasm_hic/$ID.hifiasm.hic.hap1.fa 
gfatools gfa2fa hifiasm_hic_out/*.hap2.p_ctg.gfa > asm_hifiasm_hic/$ID.hifiasm.hic.hap2.fa 
samtools faidx asm_hifiasm_hic/$ID.hifiasm.hic.hap2.fa 
seqkit stat -Ta asm_hifiasm_hic/$ID.hifiasm.hic.hap1.fa asm_hifiasm_hic/$ID.hifiasm.hic.hap2.fa > asm_hifiasm_hic/2hap.stat 2>>"$ID.hifi.log"

# Move and clean up files to organize the output
echo "Moving and cleaning up files for $ID..." >> "$ID.hifi.log"
#mv hifiasm_hic_out/*.hap*.p_ctg.gfa .
rm -rf hifiasm_hic_out
echo "Cleanup completed for $ID." >> "$ID.hifi.log"

echo "Process completed for $ID at $(date)" >> "$ID.hifi.log"