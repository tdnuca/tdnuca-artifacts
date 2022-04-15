#!/bin/bash

#rm bflog
#rm bllog
#rm cplog
size=8192
bsize=1024

#size=1024
#bsize=128
#rm bf_${size}\ ${bsize}\ 0\ file\ 1.txt
#echo "total fast avg stddev min max mean" > bf_${size}\ ${bsize}\ 0\ file\ 1.txt
#./nosmt.sh bf ${size} ${bsize} 1 #size, block_size, iterations

#echo " " > heft_${size}\ ${bsize}\ 0\ file\ 1.txt
#echo "total fast avg stddev min max mean" > heft_${size}\ ${bsize}\ 0\ file\ 1.txt
#./nosmt.sh heft ${size} ${bsize} 1 #size, block_size, iterations


rm botlev_${size}\ ${bsize}\ 0\ file\ 1.txt
echo "total fast avg stddev min max mean" > botlev_${size}\ ${bsize}\ 0\ file\ 1.txt
./nosmt.sh botlev ${size} ${bsize} 1
rm cpath_${size}\ ${bsize}\ 0\ file\ 1.txt
echo "total fast avg stddev min max mean" > cpath_${size}\ ${bsize}\ 0\ file\ 1.txt
./nosmt.sh cpath ${size} ${bsize} 1
#./pexcel.sh



