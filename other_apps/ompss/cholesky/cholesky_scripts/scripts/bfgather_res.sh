#!/bin/bash

  export NX_SCHEDULE=bf
  ./sim_chol.sh
  rm log
  ./excel.sh
  cp log bf_results.txt
  cat bf_results.txt
  rm log

