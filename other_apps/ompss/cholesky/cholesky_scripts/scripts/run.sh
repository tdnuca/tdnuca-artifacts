#!/bin/bash

matrix_size=16384

block_size=4096

#@wall_clock_limit=01:00:00

#for policy in bf dbf wf cilk priority smartpriority affinity affinity-smartpriority 
#for policy in priority smartpriority bf
#do
	for cores in 1 2 4 8 16 32
#`seq 16`
	do
		export OMP_NUM_THREADS=$cores
		#export NX_ARGS="--schedule=$policy"
		#export NX_ARGS="--schedule=botlev --update-freq=1000"
		#export NX_ARGS="--schedule=botlev"
		#export NX_ARGS="--schedule=bf --schedule-priority"
		export EXTRAE_COUNTERS=rapl:::PP0_ENERGY:PACKAGE0
                export NX_ARGS="--schedule=bf"
		#export NX_INSTRUMENTATION=extrae
		echo $NX_ARGS
		#export NX_EXTRAE_FILE_NAME=cholesky_bf_$(matrix_size)_$(block_size)_$(cores)
		#./cholesky 256 2 0 file > cholesky_bf_256_2_$cores
		./cholesky ${matrix_size} ${block_size} 0 file > cholesky_bf_${matrix_size}_${block_size}_$cores
	done
#done



        for cores in 1 2 4 8 16 32
        do
                export OMP_NUM_THREADS=$cores
                #export NX_ARGS="--schedule=$policy"
                export NX_ARGS="--schedule=botlev"
                #export NX_ARGS="--schedule=botlev"
                #export NX_ARGS="--schedule=bf --schedule-priority"
                #export NX_ARGS="--schedule=botlev"
                #export NX_INSTRUMENTATION=extrae
                echo $NX_ARGS
               # ./cholesky 256 2 0 file > cholesky_botlev_256_2_$cores
		 #export NX_EXTRAE_FILE_NAME=cholesky_bf_${matrix_size}_${block_size}_${cores}
		./cholesky  ${matrix_size} ${block_size} 0 file > cholesky_botlev_${matrix_size}_${block_size}_$cores
        done


        for cores in 1 2 4 8 16 32
        do
                export OMP_NUM_THREADS=$cores
                #export NX_ARGS="--schedule=$policy"
                #export NX_ARGS="--schedule=botlev --update-freq=1000"
                #export NX_ARGS="--schedule=botlev"
                export NX_ARGS="--schedule=bf --schedule-priority"
                #export NX_ARGS="--schedule=bf"
                #export NX_INSTRUMENTATION=extrae
                echo $NX_ARGS
               # ./cholesky 256 2 0 file > cholesky_priority_256_2_$cores
		# export NX_EXTRAE_FILE_NAME=cholesky_bf_${matrix_size}_${block_size}_${cores}
		./cholesky  ${matrix_size} ${block_size} 0 file > cholesky_priority_${matrix_size}_${block_size}_$cores
        done


