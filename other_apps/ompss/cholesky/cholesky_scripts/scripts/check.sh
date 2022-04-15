#!/bin/bash

export NX_ARGS="--schedule=botlev"
export OMP_NUM_THREADS=1
./cholesky 512 32 0 file
sort tasks.txt > dirtasks/tasks1.txt
for cores in 2 4 8 16 32
do
   export OMP_NUM_THREADS=$cores
   for s in `seq 10` 
   do
      ./cholesky 512 32 0 file 
      sort ./tasks.txt > dirtasks/tasks${cores}_${s}.txt
   done
done
#echo $(source ./diff.sh)
