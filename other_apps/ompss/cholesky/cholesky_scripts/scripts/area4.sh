#!/bin/bash

#for policy in bf
#do
#  export NX_SCHEDULE=${policy}

  export NX_HP_FROM=0
  export NX_HP_TO=3
  taskset -c 0-3 ./cholesky 8192 256 0 file 10 > res_4_0.txt

  export NX_HP_FROM=2
  export NX_HP_TO=4
  taskset -c 0-4,16-17 ./cholesky 8192 256 0 file 10 > res_3_4.txt

  export NX_HP_FROM=4
  export NX_HP_TO=5
  taskset -c 0-5,16-19 ./cholesky 8192 256 0 file 10 > res_2_8.txt

  export NX_HP_FROM=6
  export NX_HP_TO=6
  taskset -c 0-6,16-21 ./cholesky 8192 256 0 file 10 > res_1_12.txt
#done


