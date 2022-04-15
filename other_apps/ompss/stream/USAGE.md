# stream
Original code developed by John D. McCalpin
Author (StarSs port): Rosa Badia

## Issues
Current task granularity might not be optimal.

## Compilation instructions
Use any of the targets (perf, instr, debug, seq) to build the corresponding executable, e.g. `make perf` will build the default ompss executable.

## Running instructions
Pick a size for data in memory and run `./stream size`

It is the array size per node, thus with 3 arrays the total memory used is 3 * size * number of nodes.
If the size is not given, it defaults to 512M.


## Additional information
Currently the task size is always the array size split equally in 128 blocks.
If compiled with VALIDATE defined, the results of the operations performed will be checked for correctness.

**ONLY THE SMP KERNELS ARE TESTED**


```
> NX_ARGS=--smp-workers=4\ --summary ./stream $(( 256 * 1024 * 1004 ))
MSG: [?] Nanos++ Initial Environment Summary
==========================================================
===================== Global Summary =====================
=== Nanos++ version:     0.15a
=== PID:                 146416
=== Num. worker threads: 4
=== System CPUs:         active[ 2-5 ] - inactive[ 0-1, 6-47 ]
=== Binding:             true
=== Prog. Model:         OmpSs
=== Priorities:          Not needed / Disabled
=== Plugin:              SMP PE Plugin
===  | PEs:              48
===  | Worker Threads:   4
==========================================================

threads is 4, block size 2056192
-------------------------------------------------------------
STREAM version $Revision: 5.8 $
-------------------------------------------------------------
This system uses 8 bytes per DOUBLE PRECISION word.
-------------------------------------------------------------
Array size = 263192576, Offset = 0
Total memory required = 6024.0 MB.
Each test is run 10 times, but only
the *best* time for each is used.
-------------------------------------------------------------
Printing one line per active node...
node 0: a=0x2b4988000010 b=0x2b4a05801010 c=0x2b4a83002010, using regular malloc
solve_time:5230407
Function      Rate (MB/s)   Avg time     Min time     Max time
Copy:       43165.3114       0.0980       0.0976       0.0987
Scale:      43566.8190       0.0971       0.0967       0.0974
Add:        45650.5288       0.1387       0.1384       0.1390
Triad:      45564.2462       0.1395       0.1386       0.1409
-------------------------------------------------------------
TOTAL time (including initialization) =       5.2304 seconds
-------------------------------------------------------------
========== STREAM TRIAD RESULTS =========
  Execution time (sec): 0.138631
  Performance (GFLOPS): 45.564246
  Memory alocated (GB): 5.88
=========================================
-------------------------------------------------------------
MSG: [?] Nanos++ Final Execution Summary
==========================================================
=== Application ended in 5 seconds
=== 5288 tasks have been executed
==========================================================
```
