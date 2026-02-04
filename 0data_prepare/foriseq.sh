#!/bin/sh

conda activate /public/home/WHC_zhangx/anaconda3/envs/iseq 

wd=/public/home/WHC_zhangx/inversion_proj/0data/DToL_2hifiasm
cd $wd

for i in `cat ids`
do 
cd $wd/$i
mkdir hifi 
cd hifi
cat ../*.PACBIO.SRAs | while read hifiRun; do
    iseq -i ${hifiRun} -g -a &
done 

cd $wd/$i
mkdir hic 
cd hic
cat ../*.Hi-C.SRAs | while read hicRun; do
   nohup iseq -i ${hicRun} -g -a &
done 
cd $wd
done
