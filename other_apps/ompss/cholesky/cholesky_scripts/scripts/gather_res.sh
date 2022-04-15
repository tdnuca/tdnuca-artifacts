#!/bin/bash

export NX_SCHEDULE=botlev
export NX_MAXB=1
for steal in 0 1
  do
  export NX_STEALB=$steal
  for strict in 0 1
  do
    export NX_STRICTB=$strict
    ./sim_chol.sh
    rm log
    ./excel.sh
    echo "STEAL "${steal} " STRICT " ${strict}
    #cat log
    cp log botlev_steal${steal}_strict${strict}.txt
  done
done

for steal in 0 1
do
  for strict in 0 1
  do 
   echo "STEAL "${steal} " STRICT " ${strict}
   cat botlev_steal${steal}_strict${strict}.txt
  done
done
