#!/bin/bash
export NX_HP_FROM=3
export NX_HP_TO=3
export NX_MAXB=0
export NX_STRICTB=0
export NX_INSTRUMENTATION=extrae
for sched in botlev 
do
  export NX_SCHEDULE=$sched
  for cores in 4 8 16 20 24 28 30 32
  do
     for bs in 128 
     do
       for it in 1
       do
         export EXTRAE_DIR=/scratch/Computational/kchronak/cholesky_overheads/${sched}/${cores}
         NX_PES=${cores} ./cholesky_instr 1024 ${bs} 0 file ${it} 
         /apps/CEPBATOOLS/extrae/2.4.2/host/bin/mpi2prv -f /scratch/Computational/kchronak/cholesky_overheads/${sched}/${cores}/TRACE.mpits -o /scratch/Computational/kchronak/cholesky_overheads/${sched}/${cores}/cholesky_${sched}_1K_${bs}_${it}it_${cores}cores_32fr.prv
       done
     done
     rm -r /scratch/Computational/kchronak/cholesky_overheads/${sched}/${cores}/set-0
     rm -r /scratch/Computational/kchronak/cholesky_overheads/${sched}/${cores}/TRACE.mpits
  done
done

#export NX_HP_FROM=16
#export NX_HP_TO=31
#~/tasksim/installtasksim/bin/burstsim tasksimConfs/hetero3216.conf cholesky.ts > res_32_16.txt


