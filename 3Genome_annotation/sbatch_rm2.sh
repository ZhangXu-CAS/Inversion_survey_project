#!/bin/bash
#SBATCH --job-name=RM2_1026
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=120G
#SBATCH --time=60:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-17
#SBATCH --output=RM2_%A_%a.out
module load repeatmodeler/2.0.7

i=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" rm2_ids)
echo "Starting RepeatModeler pipeline at $(date) for $i.hap2"

WD=/project/6004758/zhangxu/TE_anno
INPUT=/project/6004758/zhangxu/4Annotation_done/${i}/${i}_hap2_renchr.fa
TMPDIR=$SLURM_TMPDIR/${i}
OUTDIR=${WD}/${i}
mkdir -p $TMPDIR $OUTDIR
cd $TMPDIR

cp $INPUT .

BuildDatabase -name $i.hap2 $INPUT
RepeatModeler -database $i.hap2 -engine ncbi -threads 30 &> $i.rm2.out
RepeatMasker -lib ${i}.hap2-families.fa -e ncbi -dir . ${i}_hap2_renchr.fa -pa 30 -gff

cp ${i}.hap2-families.fa ${i}_hap2_renchr.fa.out.gff ${i}_hap2_renchr.fa.masked ${i}_hap2_renchr.fa.tbl $i.rm2.out $OUTDIR/

echo "Cleaning temporary files..."
rm -rf $TMPDIR

echo "RM completed at $(date)"