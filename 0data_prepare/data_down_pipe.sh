#!/bin/bash
# Script for automatically downloading data from SRA using project IDs---Xu Zhang, Nov.04, 2024.
# Improved for readability, error handling, and efficiency by ChatGPT.
# Run at 2025-2-14
# Set working directory
wd=$(pwd)
cd "$wd" || exit 1  # Exit if directory doesn't exist

#Step 1: Create and process directories for each project ID
for id in  $(cat ids); do

    mkdir -p "$id"
    cd "$id" || continue
    
    echo "Fetching run info for $i..."
    esearch -db sra -query "$id" | efetch -format runinfo | \
    cut -d, -f1,2,5,8,13,14,19,21,22,25,28,29 > "${id}.runinfo.csv"
    
    # Filter for Hi-C and PACBIO data and save SRR accessions
    grep "Hi-C" "${id}.runinfo.csv" | cut -d, -f1 > "${id}.Hi-C.SRAs"
    grep "PACBIO" "${id}.runinfo.csv" | cut -d, -f1 > "${id}.PACBIO.SRAs"
    
    cd "$wd"
done

#Step 2: Create prefetch scripts for PACBIO and Hi-C SRAs and run them in background
for id in $(cat ids); do
    cd "$wd/$id" || continue
    
    # PACBIO prefetch script
	echo "Prefetching PACBIO data for $i..."
    cat << EOF > "${id}.hifi.prefetch.sh"
#!/bin/bash
cd $id
prefetch --option-file ${id}.PACBIO.SRAs -X 1000G
EOF
nohup bash "${id}.hifi.prefetch.sh" &

    # Hi-C prefetch script
	 echo "Prefetching Hi-C data for $i..."
    cat << EOF > "${id}.hic.prefetch.sh"
#!/bin/bash
cd $id
prefetch --option-file ${id}.Hi-C.SRAs -X 1000G
EOF
   nohup bash "${id}.hic.prefetch.sh" &
    
    cd "$wd" || exit 1
done

#Step 3: Convert SRA files to FASTQ for PACBIO
for id in $(cat ids); do
    cd "$wd/$id" || continue
	echo "Converting PACBIO SRA to FASTQ for $id..."
    mkdir -p hifi
    while read -r sra; do
        mv "$sra"/*.sra hifi/ 2>/dev/null
    done < "${id}.PACBIO.SRAs"
#	nohup fastq-dump --gzip --split-3 hifi/*.sra -O hifi &
    cd "$wd" || exit 1
done

#Step 4: Convert SRA files to FASTQ for Hi-C
for id in $(cat ids); do
	echo "Converting Hi-C SRA to FASTQ for $i..."
    cd "$wd/$id" || continue
    mkdir -p hic
	while read -r sra; do
        mv "$sra"/*.sra hic/ 2>/dev/null
    done < "${id}.Hi-C.SRAs"
#    nohup fastq-dump --gzip --split-3 hic/*.sra -O hic &
    cd "$wd" || exit 1
done

#Step 5: Concatenate mutiple Hi-C FASTQ files
for mid in $(cat mlti.ids); do
    {
        cd "$wd/$mid/hic" || continue
        cat *_1.fastq.gz > "${mid}.hic_1.fastq.gz"
        cat *_2.fastq.gz > "${mid}.hic_2.fastq.gz"
        cd "$wd" || exit 1
    } &
done

wait  # Wait for all background jobs to finish before exiting script
echo "Download and processing completed."
