#!/bin/bash


for steal in 0 1
do

  for strict in 0 1
  do
    export NX_SCHEDULE=botlev
    export NX_STEALB=${steal}
    export NX_STRICTB=${strict}
    echo "--- botlev no strict simple stealing ---"
    echo "--- area4: ---" 
    ./area4.sh
    ./pexcel.sh
    echo "--- area8: ---"
    rm res_*
    ./area8.sh
    ./pexcel.sh
  done
done

cp -r /home/Computational/kchronak/ompss/bf_install/* /home/Computational/kchronak/ompss/installnanos/

export NX_SCHEDULE=bf
echo "--- bf ---"
echo "--- area4: ---"
rm res_*
./area4.sh
./pexcel.sh
echo "--- area8: ---"
rm res_*
./area8.sh
./pexcel.sh



