#!/bin/bash
#SBATCH --job-name=1220bamyahs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-67

# General configuration
threads=24
wd="/lustre07/scratch/zhangxu/inversion_proj/1220_bamyahs"
module load java/21.0.1

# Read the ID for this task from the file 'ids'
id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ids)
id_dir="$wd/${id}"
log_file="$id_dir/${id}_yahs.log"

# Validate ID and working directory
if [[ -z $id || ! -d "$id_dir" ]]; then
    echo "Error: ID is empty or working directory does not exist for $id" >> "$log_file"
    exit 1
fi

pushd "$id_dir" > /dev/null || exit 1
echo "$id - $(date)" | tee -a "$log_file"

# Step 1: Index reference assemblies
for hap in hap1 hap2; do
    ref_file="asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa"
    bwt_file="${ref_file}.bwt"
    if [[ ! -f "$bwt_file" || ! -s "$bwt_file" ]]; then
        echo "Indexing reference for $hap..." | tee -a "$log_file"
        bwa index "$ref_file" > "$log_file" 2>&1 &
    else
        echo "Reference already indexed for $hap, skipping." > "$log_file"
    fi
done
wait

# Step 2: Align Hi-C reads and filter BAM files
bam_status=$(samtools quickcheck ${id}.hap*.HiC.filtered.bam 2>/dev/null && echo "valid" || echo "invalid")
if [[ $bam_status != "valid" ]]; then
    for hap in hap1 hap2; do
        ref_file="asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa"
        bam_file="${id}.${hap}.HiC.bam"
        filtered_bam="${id}.${hap}.HiC.filtered.bam"

        if [[ ! -f "$bam_file" || ! -s "$bam_file" ]]; then
            echo "Aligning reads to generate BAM file for $hap..." >> "$log_file"
            bwa mem -5SP -t "$threads" "$ref_file" hic_data/*_1.fastq.gz hic_data/*_2.fastq.gz | \
            samblaster | samtools view -@ "$threads" -b -F 3340 -o "$bam_file" >> "$log_file" 2>&1
        else
            echo "BAM file already exists for $hap, skipping alignment." >> "$log_file"
        fi

        if [[ ! -f "$filtered_bam" || ! -s "$filtered_bam" ]]; then
            echo "Filtering BAM file for $hap..." >> "$log_file"
            filter_bam "$bam_file" 1 --nm 3 --threads "$threads" | \
            samtools view -b -@ "$threads" -o "$filtered_bam" >> "$log_file" 2>&1
        else
            echo "Filtered BAM file already exists for $hap, skipping filtering." >> "$log_file"
        fi

        if [[ -s "$filtered_bam" ]]; then
            rm "$bam_file" && echo "Removed HiC.bam file to save disk space" >> "$log_file"
        fi
    done
fi

# Step 3: Run Yahs scaffolding
for hap in hap1 hap2; do
    output_dir="${hap}_yahs"
    mkdir -p "$output_dir"
    echo "Processing Yahs for $hap..." | tee -a "$log_file"
    yahs "asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa" "${id}.${hap}.HiC.filtered.bam" -o "$output_dir/${id}.${hap}" >> "$log_file" 2>&1 &
done
wait
echo "Yahs processing complete for $id - $(date)" | tee -a "$log_file"

# Step 4: Extract scaffolding results
juicer_tool="/home/zhangxu/software/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar"
for hap in hap1 hap2; do
    chr_dir="${hap}_yahs"
    if [[ -d "$chr_dir" ]]; then
        (
            cd "$chr_dir" || { echo "Failed to enter $chr_dir"; exit 1; }

            echo "Processing $hap started at $(date)" >> "$log_file"

            # Generate Juicer TXT
            if juicer pre -a -o ${id}_${hap} \
                "$id_dir/${id}.${hap}.HiC.filtered.bam" \
                "${id}.${hap}_scaffolds_final.agp" \
                "$id_dir/asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa.fai" > "${id}_${hap}_juicer.log" 2>&1; then
                echo "Juicer TXT generation complete for $hap" >> "$log_file"
            else
                echo "Juicer TXT generation failed for $hap" >> "$log_file"
                continue
            fi

            # Generate Juicer HIC
            if java -jar "${juicer_tool}" \
                pre ${id}_${hap}.txt ${id}_${hap}.hic.part \
                <(grep PRE_C_SIZE ${id}_${hap}_juicer.log | awk '{print $2" "$3}') >> "$log_file" 2>&1; then
                mv "${id}_${hap}.hic.part" "${id}_${hap}.hic"
                echo "Generated .hic file for $hap: ${id}_${hap}.hic" >> "$log_file"
            else
                echo "Failed to generate .hic file for $hap" >> "$log_file"
            fi
        ) &
    else
        echo "Skipping $hap as $chr_dir does not exist" >> "$log_file"
    fi
done
wait
echo "All tasks complete for $id - $(date)" | tee -a "$log_file"
popd > /dev/null
