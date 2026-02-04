#!/bin/bash
## For BUSCO evaluation of assembled genomes (plants, animals, fungi)
## Xu Zhang -- 2025-08-14

conda activate busco
WD=/Volumes/70T/AsteroidScratch/xuzhang/inver_proj/5BUSCO

run_busco() {
    local kingdom_dir=$1
    local lineage=$2

    cd "$WD/$kingdom_dir" || { echo "Directory $WD/$kingdom_dir not found"; exit 1; }
    for i in $(cat ids); do
        busco -i "${i}_hap1_revchr.fa" \
              -l "$lineage" \
              -m geno \
              -o "$i" \
              -c 20 \
              --offline
    done
}

# Run the three kingdoms in parallel
run_busco "Animal" "./metazoa_odb12" &     # 动物
run_busco "Plants" "./viridiplantae_odb12" & # 植物
run_busco "Funga"  "./fungi_odb12" &         # 真菌

wait
echo "All BUSCO runs completed at $(date)"

for i in `cat ids`; do
python parse_busco_results_single.py $i $i/*.txt
done