for i in `cat ids`;
do 
python get_inv_info.py $i
echo "get_inv_info done for $i at `date`"
done
