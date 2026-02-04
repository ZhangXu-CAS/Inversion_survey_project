##For simulation to detect accuracy of minimap + syri pipeline 
## We have 23 genome, 10 from animal, 10 from plants and 3 from fungal
## Using three length and counts parameters: N=80, 60, 40. L< 10Kb, 100Kb, 1Mb
## to reprsent: many small inversion, few big inversion and middle.
for i in `cat ids`; do
mkdir -p Re_ani/$i Re_mm/$i Re_pla/$i
SURVIVOR simSV  Allgeno/${i}_hap1_revchr.fa par_ani.txt 0.0001 0 Re_ani/$i/$i.sim
SURVIVOR simSV  Allgeno/${i}_hap1_revchr.fa par_mid.txt 0.0001 0 Re_mm/$i/$i.sim
SURVIVOR simSV  Allgeno/${i}_hap1_revchr.fa par_plants.txt 0.0001 0 Re_pla/$i/$i.sim
echo "$i sim done at `date`"
done

cd /Volumes/70T/AsteroidScratch/xuzhang/inver_proj/7Sim/Re_ani
for i in `cat ids`; do
SURVIVOR eval $i/syri_results/syri.vcf  $i/$i.sim.bed 10 $i.eval_res
done
