#!/bin/bash

for cores in 1 2 4 6 8 12 16 20 24 28 30 32
do
  for big in 1 2 4 6 8 10 12 16
  do
#    eval from=`echo "$cores - $big" | bc -l`
#    eval to=`echo "$cores - 1" | bc -l`
#    echo ${from} " - " ${to} " for " ${cores} " cores and " ${big} " HP "
#    export NX_HP_FROM=$from
#    export NX_HP_TO=$to

    if [ ${cores} -gt ${big} ] 
      then 
        eval from=`echo "$cores - $big" | bc -l`
        eval to=`echo "$cores - 1" | bc -l`
        echo ${from} " - " ${to} " for " ${cores} " cores and " ${big} " HP "
        export NX_HP_FROM=$from
        export NX_HP_TO=$to

        ~/tasksim/installtasksim/bin/burstsim tasksimConfs/hetero${cores}${big}.conf cholesky.ts > res_${cores}_${big}.txt 
    fi
    if [ ${cores} -eq 1 -a ${big} -eq 1 ] 
      then 
        eval from=`echo "$cores - $big" | bc -l`
        eval to=`echo "$cores - 1" | bc -l`
        echo ${from} " - " ${to} " for " ${cores} " cores and " ${big} " HP "
        export NX_HP_FROM=$from
        export NX_HP_TO=$to
        ~/tasksim/installtasksim/bin/burstsim tasksimConfs/hetero${cores}${big}.conf cholesky.ts > res_${cores}_${big}.txt 
    fi
  done
done

#export NX_HP_FROM=16
#export NX_HP_TO=31
#~/tasksim/installtasksim/bin/burstsim tasksimConfs/hetero3216.conf cholesky.ts > res_32_16.txt


