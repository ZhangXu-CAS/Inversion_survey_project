#!/bin/bash
#SBATCH --job-name=RM1_1026
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=120G
#SBATCH --time=60:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-18
#SBATCH --output=RM1_%A_%a.out
module load repeatmodeler/2.0.7

i=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" rm_ids)
echo "Starting RepeatModeler pipeline at $(date) for $i.hap1"

# 设置路径
WD=/project/6004758/zhangxu/TE_anno
INPUT=/project/6004758/zhangxu/4Annotation_done/${i}/${i}_hap1_revchr.fa
TMPDIR=$SLURM_TMPDIR/${i}
OUTDIR=${WD}/${i}
mkdir -p $TMPDIR $OUTDIR
cd $TMPDIR

cp $INPUT .

BuildDatabase -name $i.hap1 $INPUT
RepeatModeler -database $i.hap1 -engine ncbi -threads 30 &> $i.rm1.out
RepeatMasker -lib ${i}.hap1-families.fa -e ncbi -dir . ${i}_hap1_revchr.fa -pa 30 -gff

cp ${i}.hap1-families.fa ${i}_hap1_revchr.fa.out.gff ${i}_hap1_revchr.fa.masked ${i}_hap1_revchr.fa.tbl $i.rm1.out $OUTDIR/

echo "Cleaning temporary files..."
rm -rf $TMPDIR

echo "RM completed at $(date)"