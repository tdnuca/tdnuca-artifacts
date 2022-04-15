#!/bin/bash


      grep "time (s)"  bf_*_*.txt | awk -F":" '{print $1, $3}'  >>bflog;

sed -i 's/bf_//g' bflog; sed -i 's/\.txt//g' bflog; sed -i 's/_/ /g' bflog

sort -r timelog


      grep "performance (gflops)"  res_*_*.txt | awk -F":" '{print $1, $3}'  >>perflog;
sed -i 's/res_//g' perflog; sed -i 's/\.txt//g' perflog; sed -i 's/_/ /g' perflog
sort -r perflog

