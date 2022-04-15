# specfem3D
Specfem3D Globe Version 4.0 - A high-order finite-element earthquake modeling application  
Author: Dimitri Komatitsch and Jeroen Tromp, ported to OmpSs by Unknown (Imported from HCA-Apps repository)  
Contact: Vladimir Dimic <vladimir.dimic@bsc.es>  

## Issues
With module ompss/stable (PATH): OmpSs 16.06.3 - nanox 0.10.3, mcxx 2.0.0, following warning appears:  
`specfem3D_smpss_v0.81.c: In function ‘main’:`  
`specfem3D_smpss_v0.81.c:3270:35: warning: initialization from incompatible pointer type`  
`specfem3D_smpss_v0.81.c:3281:15: warning: passing argument 3 of ‘record_seis’ from incompatible pointer type`  
`specfem3D_smpss_v0.81.c:671:6: note: expected ‘float (*)[3]’ but argument is of type ‘float *’`

## Compilation instructions
Run "make" to generate native binary  
Run "make instr" to generate instrumented binary (useful for e.g. tracing for TaskSim simulations)  

## Running instructions
`./specfem3D`  
The application does not accept any command line arguments.  

## Inputs
Two input files are necessary for running the application. They need to be in the same folder as the binary. The files can be downloaded from:  
`/home/bsc18/bsc18902/scratch/soft/test/hca-apps/bin/matrices.dat`  
`/home/bsc18/bsc18902/scratch/soft/test/hca-apps/bin/proc000000_reg1_database.dat`  

In the file `values_from_mesher_C.h`, there are numerical constants that can be changed.  

## Additional information
###  Environment info
Vladimir (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7)  
TaskSim tracing done with (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7). For details ask Vladimir.  

### Output with --summary
`MSG: [?] ========== Nanos++ Initial Environment Summary ==========`  
`MSG: [?] === PID:                 22692`  
`MSG: [?] === Num. worker threads: 16`  
`MSG: [?] === System CPUs:         active[ 0-15 ] - inactive[  ]`  
`MSG: [?] === Binding:             true`  
`MSG: [?] === Prog. Model:         OmpSs`  
`MSG: [?] === Priorities:          Not needed / disabled`  
`MSG: [?] === Plugin:              SMP PE Plugin`  
`MSG: [?] ===  | PEs:              16`  
`MSG: [?] ===  | Worker Threads:   16`  
`MSG: [?] === Plugin:              OpenCL PE Plugin`  
`MSG: [?] ===  | PEs:              0`  
`MSG: [?] ===  | Worker Threads:   0`  
`MSG: [?] =========================================================`  
`MSG: [?] ============ Nanos++ Final Execution Summary ============`  
`MSG: [?] === Application ended in 6 seconds`  
`MSG: [?] === 2130 tasks have been executed`  
`MSG: [?] =========================================================`  

