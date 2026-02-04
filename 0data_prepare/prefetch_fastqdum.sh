#!/bin/bash

wd=/home/zhangxu/inversion_proj/0data/DToL_SRAs
cd $wd

IDs=$(cat bioprojIDs)
for i in $IDs
do
cd $wd/$i
esearch -db sra -query $i | efetch -format runinfo | cut -d, -f1,2,5,8,13,14,19,21,22,25,28,29 > $i.runinfo.csv
cat $i.runinfo.csv | grep Hi-C | cut -d, -f1 > $i.Hi-C.SRAs
cat $i.runinfo.csv | grep PACBIO | cut -d, -f1 > $i.PACBIO.SRAs
prefetch --option-file $i.Hi-C.SRAs -X 500G &
prefetch --option-file $i.PACBIO.SRAs -X 500G &
pacsra=$(cat ${i}.PACBIO.SRAs)
for j in $pacsra
do

mv $j $j.hifi 
cd $j.hifi
fastq-dump --gzip --split-3 *.sra &
cd $wd/$i
done

hiccsra=$(cat ${i}.Hi-C.SRAs)
for k in $hiccsra
do
mv $k $k.hic 
cd $k.hic
fastq-dump --gzip --split-3 *.sra 
cd $wd/$i
done
cd $wd
done
