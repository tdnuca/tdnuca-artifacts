#!/bin/bash


#grep "Final cycle count"  res_28*.txt | awk -F":" '{print $1, $3}'  >log_28; sed -i 's/res_//g' log_28; sed -i 's/\.txt//g' log_28; sed -i 's/_/ /g' log_28 >>log
for cores in 1 2 4 6 8 12 16 20 24 28 30 32
do
#  for big in 1 2 4 6 8 16
#  do 
#    if [ ${cores} -gt ${big} ] 
#      then
      grep "Final cycle count"  res_${cores}_*.txt | awk -F":" '{print $1, $3}'  >>log;
#    fi
#  done
done

sed -i 's/res_//g' log; sed -i 's/\.txt//g' log; sed -i 's/_/ /g' log 







