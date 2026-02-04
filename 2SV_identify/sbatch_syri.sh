#!/bin/bash
#SBATCH --job-name=0929mapping
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-1

module load python/3.13.2
source  /home/zhangxu/anaconda3/bin/activate base

cd /lustre07/scratch/zhangxu/inversion_proj/plant_published
# Read the ID for this task from the file 'ids'
id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ids)

bash minimap_syri.sh $id
echo "$i completed at `date`"
done

