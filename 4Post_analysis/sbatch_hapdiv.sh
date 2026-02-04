#!/bin/bash
#SBATCH --job-name=1006Hapdiv
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-6
#SBATCH --out=sbatch_Hapdiv.log
module load minimap2/2.28

wd="/project/6004758/zhangxu/Hap_div"
cd $wd
# Read sample ID from list
id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ids2)
id_dir="$wd/${id}"

cd "$id_dir" 
echo "Calculate hap_divergence for $i at $(date)"

#minimap2 -t 30 -ax asm5 --eqx ${id}_hap1* ${id}_hap2* > 2nd_chr_alignment.sam 
samtools view -bS 2nd_chr_alignment.sam | samtools sort -o $id.hap2_vs_hap1.sorted.bam
samtools index $i.hap2_vs_hap1.sorted.bam
python3 $wd/hap_divergence.py -b $id.hap2_vs_hap1.sorted.bam -o $id.hapdiv --bin 10000 --min_mapq 20

echo "Hap_divergence done for $i at $(date)"