#!/bin/bash


  export NX_HP_FROM=0
  export NX_HP_TO=8
  taskset -c 0-8 ./cholesky 8192 256 0 file 10 > res_8_0.txt

  export NX_HP_FROM=4
  export NX_HP_TO=9
  taskset -c 0-9,16-19 ./cholesky 8192 256 0 file 10 > res_6_8.txt

  export NX_HP_FROM=8
  export NX_HP_TO=11
  taskset -c 0-11,16-23 ./cholesky 8192 256 0 file 10 > res_4_16.txt

  export NX_HP_FROM=11
  export NX_HP_TO=12
  taskset -c 0-12,16-26 ./cholesky 8192 256 0 file 10 > res_2_24.txt


