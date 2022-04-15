df=`diff -r /home/Computational/kchronak/ompss/installnanos/ /home/Computational/kchronak/ompss/bf_install/`
if [[ ${df}  ]]; then
  cp -r /home/Computational/kchronak/ompss/bf_install/* /home/Computational/kchronak/ompss/installnanos/
fi
rm cholesky.ts.*
~/tasksim/installtasksim/bin/ompss ./cholesky 16384 1024 0 file 1

