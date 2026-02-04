#!/bin/bash
cd /project/6004758/zhangxu/Veri_INV

for i in `cat ids`;
do
python verify_inversions.py \
  -b ./BAM/$i/$i.sorted.bam \
  -i syri4veri/$i/$i.inv.csv -o syri4veri/$i/$i.inv.verify.csv 
done

