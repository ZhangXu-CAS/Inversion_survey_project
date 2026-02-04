#!/bin/bash

wd="/Volumes/70T/AsteroidScratch/xuzhang/inver_proj/2map"
cd "$wd"

log_file="$wd/batch_extract_chr.log"
echo "Batch processing started at $(date)" > "$log_file"

while read -r species_id; do  
    cd $species_id 
	python "$wd/extract_chr.py" "$species_id"  >> $log_file 2>&1   
    echo "Finished processing $species_id at $(date)" >> "$log_file"
cd $wd
done < ids

echo "Batch processing completed at $(date)" >> "$log_file"
