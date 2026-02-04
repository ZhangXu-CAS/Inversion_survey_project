#!/bin/bash
#SBATCH --job-name=maphifiWH_0908
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --account=rrg-rieseber-ac
#SBATCH --array=0-237
#SBATCH --out=Veri_INV_wuhan_BAMs.log

module load minimap2
LOGFILE=/project/6004758/zhangxu/Veri_INV_wuhan/Veri_INV_wuhan_BAMs.log
ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ids)
echo "[$(date)] Start processing $ID (Task $SLURM_ARRAY_TASK_ID)" >> $LOGFILE

BASE_DIR=/project/6004758/zhangxu/Veri_INV
REF=$BASE_DIR/genomes/${ID}_hap1_revchr.fa
cd $BASE_DIR

#ref index
if [ ! -f ${REF}.mmi ]; then
    echo "[$(date)] [$ID] Building minimap2 index ..." >> $LOGFILE
    minimap2 -d ${REF}.mmi $REF
fi

if [ ! -f ${REF}.fai ]; then
    echo "[$(date)] [$ID] Building samtools faidx ..." >> $LOGFILE
    samtools faidx $REF
fi

SAMPLEDIR=$BASE_DIR/BAM/$ID
mkdir -p $SAMPLEDIR

FASTQ_DIR=$BASE_DIR/hifi_reads/$ID
FASTQ_FILES=$(ls $FASTQ_DIR/*.fastq 2>/dev/null)

if [ -z "$FASTQ_FILES" ]; then
    echo "[$(date)] [ERROR][$ID] No FASTQ found in $FASTQ_DIR" >> $LOGFILE
    exit 1
fi

#BAM mapping
echo "[$(date)] [$ID] Mapping start ..." >> $LOGFILE
minimap2 -ax map-hifi -t $SLURM_CPUS_PER_TASK $REF $FASTQ_FILES \
    | samtools sort -@ $SLURM_CPUS_PER_TASK -o $SAMPLEDIR/${ID}.sorted.bam

samtools index $SAMPLEDIR/${ID}.sorted.bam

#BAM checking
if samtools quickcheck $SAMPLEDIR/${ID}.sorted.bam; then
    STATS=$(samtools flagstat -@ $SLURM_CPUS_PER_TASK $SAMPLEDIR/${ID}.sorted.bam | head -n 5 | tr '\n' ';')
    echo "[$(date)] [$ID] Mapping finished successfully. Stats: $STATS" >> $LOGFILE
else
    echo "[$(date)] [ERROR][$ID] BAM validation failed for $SAMPLEDIR/${ID}.sorted.bam" >> $LOGFILE
    exit 1
fi

echo "[$(date)] [$ID] Completed." >> $LOGFILE
