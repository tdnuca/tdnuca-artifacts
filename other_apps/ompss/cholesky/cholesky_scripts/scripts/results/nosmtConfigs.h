#!/bin/bash

#set all as slow
for core in 0 1 2 3 4 5 6 7 
do
   cpufreq-set -c $core -f 1600000
done
size=$2
block_size=$3
iter=$4
<<COMMENT1
COMMENT1
#============================================================================

for steal in 0 1
do
  for strict in 0 1 
  do
#if [[ $1="botlev" ]]; then
export NX_SCHEDULE=$1
#total 8 cores
total=8
fast=1
cpufreq-set -c 7 -f 2668000
export NX_HP_FROM=7
export NX_HP_TO=7
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-7 $1 $total $fast $steal $strict 
#taskset -c 0-7 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt
#if [[ $1="cpath" ]]; then
fast=2
cpufreq-set -c 6 -f 2668000
export NX_HP_FROM=6
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-7 $1 $total $fast $steal $strict
#taskset -c 0-7 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

fast=3
cpufreq-set -c 5 -f 2668000
export NX_HP_FROM=5
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-7 $1 $total $fast $steal $strict
#taskset -c 0-7 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

fast=4
cpufreq-set -c 4 -f 2668000
export NX_HP_FROM=4
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-7 $1 $total $fast $steal $strict
#taskset -c 0-7 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

#total 6 cores
total=6
fast=1
cpufreq-set -c 4 -f 1600000
export NX_HP_FROM=5
export NX_HP_TO=5
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-5 $1 $total $fast $steal $strict
#taskset -c 0-5 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

fast=2
cpufreq-set -c 4 -f 2668000
export NX_HP_FROM=4
export NX_HP_TO=5
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-5 $1 $total $fast $steal $strict
#taskset -c 0-5 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

fast=3
cpufreq-set -c 3 -f 2668000
export NX_HP_FROM=3
export NX_HP_TO=5
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-5 $1 $total $fast $steal $strict
#taskset -c 0-5 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

fast=4
cpufreq-set -c 2 -f 2668000
export NX_HP_FROM=2
export NX_HP_TO=5
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-5 $1 $total $fast $steal $strict
#taskset -c 0-5 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

#total 4 cores
total=4
fast=1
cpufreq-set -c 2 -f 1600000
export NX_HP_FROM=3
export NX_HP_TO=3
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-3 $1 $total $fast $steal $strict
#taskset -c 0-3 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt
#fi
fast=2
cpufreq-set -c 2 -f 2668000
export NX_HP_FROM=2
export NX_HP_TO=3
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-3 $1 $total $fast $steal $strict
#taskset -c 0-3 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

fast=3
cpufreq-set -c 1 -f 2668000
export NX_HP_FROM=1
export NX_HP_TO=3
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-3 $1 $total $fast $steal $strict
#taskset -c 0-3 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

#total 2 cores
total=2
fast=1
export NX_HP_FROM=1
export NX_HP_TO=1
./runSingle.sh "${size} ${block_size} 0 file ${iter}" 0-1 $1 $total $fast $steal $strict
#taskset -c 0-1 ../cholesky ${size} ${block_size} 0 file ${iter} > $1_${total}_${fast}.txt

done
done

#fi
