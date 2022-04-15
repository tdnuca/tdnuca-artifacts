# tinyjpeg
This application benchmark performs decoding of JPEG images with fixed encoding of 2x2 MCU size and YUV color, producing RGB .tga output files. Images can contain up to one RST marker per MCU line.  
Author: Luc Saillard <luc@saillard.org>, ported to OmpSs by Unknown (Imported from HCA-Apps repository)  
Contact: Vladimir Dimic <vladimir.dimic@bsc.es>  

## Issues
With the module OmpSs 16.06.3 - nanox 0.10.3, mcxx 2.0.0, following error occurs:  
`terminate called after throwing an instance of 'nanos::FatalError'`  
`  what():  FATAL ERROR: [1] invalid key, can not continue. Address 0 w/len 4 [wd desc. not available] conflicts with address: 0x1, size: 4`  

## Compilation instructions
Run "make" to generate native binary  
Run "make instr" to generate instrumented binary (useful for e.g. tracing for TaskSim simulations)  

## Running instructions
`./tinyjpeg [option] <input_filename.jpeg> <output_filename>`  
`option:`  
`  --benchmark - Measure the timing of 1 conversion`  

## Inputs
Input file is a .jpeg picture. An example picture can be found at  
`/home/bsc18/bsc18902/soft/opt/hca-apps/hca-apps/src/kernels/tinyjpeg/inputs/scc_reference.jpg`  

## Additional information
###  Environment info
Vladimir (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7)  
TaskSim tracing done with (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7). For details ask Vladimir.  

### Output with --summary
`MSG: [?] ========== Nanos++ Initial Environment Summary ==========`  
`MSG: [?] === PID:            10539`  
`MSG: [?] === Num. threads:   16`  
`MSG: [?] === Active CPUs:    [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, ]`  
`MSG: [?] === Binding:        true`  
`MSG: [?] === Prog. Model:    OmpSs`  
`MSG: [?] === Plugin:         OpenCL PE Plugin`  
`MSG: [?] ===  | Threads:     0`  
`MSG: [?] =========================================================`  
`MSG: [?] ============ Nanos++ Final Execution Summary ============`  
`MSG: [?] === Application ended in 1 seconds`  
`MSG: [?] === 4000 tasks have been executed`  
`MSG: [?] =========================================================`  

