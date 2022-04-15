# stap
This Space-Time Adaptive Processing (STAP) application is a simplified view of an MTI (Moving Target Indication) application, whose goal is to receive the echo from the ground of a periodic sequence of radar pulses, and to detect the objects that are moving on the ground among all the other generally still reflecting surfaces under the radar beam (ground clutter). Basically, radar processing permits to estimate both the position of a target through the delay between transmission of a signal (pulse) and reception of its echo, and its speed through the Doppler effect that affects echoes of several identical pulses sent periodically: the speed of the target results in a (small) variation of its distance from the radar, which is only visible as a phase shift on the radar signal (e.g., around 10 GHz). In this basic approach, Doppler processing consists in a bank of filters (e.g., a Fast Fourier Transform) each one tuned towards a particular phase shift between successive echoes. This application is written in C.   
Author: T. Petrisor <claudia-teodora.petrisor@thalesgroup.com>, ported to OmpSs by Unknown (Imported from HCA-Apps repository)  
Contact: Vladimir Dimic <vladimir.dimic@bsc.es>  

## Issues
With the module OmpSs 16.06.3 - nanox 0.10.3, mcxx 2.0.0, running application fails:  
`terminate called after throwing an instance of 'nanos::FatalError'`  
`what():  FATAL ERROR: [2] Error: cd.getNumDimensions() returns 2 but I already have the object registered with 3 dimensions. WD is : n/a copy index: 0 got reg object? 0`  

## Compilation instructions
Run "make" to generate native binary  
Run "make instr" to generate instrumented binary (useful for e.g. tracing for TaskSim simulations)  

## Running instructions
`./stap`  
The application does not accept any command line arguments. To be able to run the application, you need to place the file InputStimuli.txt in the same folder as the binary.  

## Inputs
Input file can be downloaded from  
`/gpfs/scratch/bsc18/bsc18902/soft/test/hca-apps/bin/InputStimuli.txt`  

## Additional information
###  Environment info
Vladimir (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7)  
TaskSim tracing done with (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7). For details ask Vladimir.  

### Output with --summary
`MSG: [?] ========== Nanos++ Initial Environment Summary ==========`  
`MSG: [?] === PID:            25059`  
`MSG: [?] === Num. threads:   16`  
`MSG: [?] === Active CPUs:    [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, ]`  
`MSG: [?] === Binding:        true`  
`MSG: [?] === Prog. Model:    OmpSs`  
`MSG: [?] === Plugin:         OpenCL PE Plugin`  
`MSG: [?] ===  | Threads:     0`  
`MSG: [?] =========================================================`  
`MSG: [?] ============ Nanos++ Final Execution Summary ============`  
`MSG: [?] === Application ended in 3 seconds`  
`MSG: [?] === 241515 tasks have been executed`  
`MSG: [?] =========================================================` 

