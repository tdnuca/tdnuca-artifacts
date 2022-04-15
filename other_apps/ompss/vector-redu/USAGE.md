# vector_redu
Computes the summation of the elements of a vector.  
Author: ported to OmpSs by Unknown (Imported from HCA-Apps repository)  
Contact: Vladimir Dimic <vladimir.dimic@bsc.es>  

## Issues

## Compilation instructions
Run "make" to generate native binary  
Run "make instr" to generate instrumented binary (useful for e.g. tracing for TaskSim simulations)  

## Running instructions
`./vector_redu <number_of_tasks>`  
The input parameter `number_of_tasks` specifies the number of tasks to create.  

## Inputs
Inside the source file, `NUM_PER_TASK` specifies the number of data to process per task.  

## Additional information
###  Environment info
Vladimir (MareNostrum): gcc (4.9.1) nanox (0.10.3) mcc (2.0.0)  
TaskSim tracing done with (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7). For details ask Vladimir.  

### Output with --summary
`MSG: [?] ========== Nanos++ Initial Environment Summary ==========`  
`MSG: [?] === PID:                 24785`  
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
`MSG: [?] === Application ended in 1 seconds`  
`MSG: [?] === 8192 tasks have been executed`  
`MSG: [?] =========================================================`  

