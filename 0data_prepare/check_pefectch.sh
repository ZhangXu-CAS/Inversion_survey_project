#!/bin/bash
# Check the integrity of downloaded SRA files for each ID.
# Re-run the prefetch command if files are incomplete.

# Set working directory
wd="/public/home/WHC_zhangx/inversion_proj/0data/DoL_add"
cd "$wd" || exit 1

# Loop through each ID from the ids file
for id in $(cat ids); do
    # Navigate to the folder for the current ID
    cd "$wd/$id" || continue
    
    # PACBIO prefetch script creation
    echo "Checking PACBIO data for $id..."
    cat << EOF > "${id}.hifi.prefetch.sh"
#!/bin/bash
cd "$wd/$id"
prefetch --option-file ${id}.PACBIO.SRAs -X 1000G
EOF
    
    # Check if PACBIO SRA files are already downloaded and complete
    pacbio_missing=false
    while read -r sra_file; do
        if [[ ! -f "$sra_file.sra" || $(stat -c%s "$sra_file.sra") -le 0 ]]; then
            pacbio_missing=true
            echo "PACBIO SRA file $sra_file is missing or incomplete for $id."
        fi
    done < "${id}.PACBIO.SRAs"
    
    # Run PACBIO prefetch script if files are incomplete
    if [[ $pacbio_missing == true ]]; then
        echo "Running PACBIO prefetch for $id..."
        nohup bash "${id}.hifi.prefetch.sh" &
    else
        echo "PACBIO data for $id is complete."
    fi

    # Hi-C prefetch script creation
    echo "Checking Hi-C data for $id..."
    cat << EOF > "${id}.hic.prefetch.sh"
#!/bin/bash
cd "$wd/$id"
prefetch --option-file ${id}.Hi-C.SRAs -X 1000G
EOF

    # Check if Hi-C SRA files are already downloaded and complete
    hic_missing=false
    while read -r sra_file; do
        if [[ ! -f "$sra_file.sra" || $(stat -c%s "$sra_file.sra") -le 0 ]]; then
            hic_missing=true
            echo "Hi-C SRA file $sra_file is missing or incomplete for $id."
        fi
    done < "${id}.Hi-C.SRAs"

    # Run Hi-C prefetch script if files are incomplete
    if [[ $hic_missing == true ]]; then
        echo "Running Hi-C prefetch for $id..."
        nohup bash "${id}.hic.prefetch.sh" &
    else
        echo "Hi-C data for $id is complete."
    fi

    # Return to the working directory
    cd "$wd" || exit 1
done
