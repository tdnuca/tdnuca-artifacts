# rgbyuv
Rotates image by a specified angle and converts color space from RGB to YUV
Author: Michael Andersch (Imported from HCA-Apps repository)
Contact: Vladimir Dimic (vladimir.dimic@bsc.es)

## Issues

## Compilation instructions
Run "make" to generate native binary
Run "make instr" to generate instrumented binary (useful for e.g. tracing for TaskSim simulations)

## Running instructions
./rgbyuv <input_image> <output_image> <angle> <num_threads>
<input_image>, <output_image>: source and destination image in .ppm format
<angle>: how much to rotate the input picture (in degrees)
<num_threads>: number of threads

## Inputs
Input can be any picture in .ppm format.
Example: /home/bsc18/bsc18902/soft/opt/hca-apps/hca-apps/src/kernels/tinyjpeg/inputs/llama.ppm

## Additional information
###  Environment info
Vladimir (MareNostrum): gcc (4.9.1) nanox (0.10.13) mcc (2.0.0)
TaskSim tracing done with (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7). For details ask Vladimir.

### Output with --summary

MSG: [?] ========== Nanos++ Initial Environment Summary ==========
MSG: [?] === PID:                 8572
MSG: [?] === Num. worker threads: 16
MSG: [?] === System CPUs:         active[ 0-15 ] - inactive[  ]
MSG: [?] === Binding:             true
MSG: [?] === Prog. Model:         OmpSs
MSG: [?] === Priorities:          Not needed / disabled
MSG: [?] === Plugin:              SMP PE Plugin
MSG: [?] ===  | PEs:              16
MSG: [?] ===  | Worker Threads:   16
MSG: [?] === Plugin:              OpenCL PE Plugin
MSG: [?] ===  | PEs:              0
MSG: [?] ===  | Worker Threads:   0
MSG: [?] =========================================================
MSG: [?] ============ Nanos++ Final Execution Summary ============
MSG: [?] === Application ended in 3 seconds
MSG: [?] === 9898 tasks have been executed
MSG: [?] =========================================================

