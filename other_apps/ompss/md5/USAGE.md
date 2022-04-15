# md5
This is an OpenSSL-compatible implementation of the RSA Data Security, Inc.
MD5 Message-Digest Algorithm (RFC 1321).
Check the .c files' headers for more information.
Author: Alexander Peslyak, ported to OmpSs by Michael Andersch (Imported from HCA-Apps repository)
Contact: Vladimir Dimic (vladimir.dimic@bsc.es)

## Issues
With current OMPSS module in MareNostrum (load ompss/stable (PATH): OmpSs 16.06.3 - nanox 0.10.3, mcxx 2.0.0), the application throws Nanos::FatalError
what():  FATAL ERROR: [0] invalid key, can not continue. Address 0x6bc1b0 w/len 16384 [wd desc. not available] conflicts with address: 0, size: 0

## Compilation instructions
Run "make" to generate native binary
Run "make instr" to generate instrumented binary (useful for e.g. tracing for TaskSim simulations)

## Running instructions
Usage: ./md5 <options>
-i inputnum          choose input data set
-c iterations        specify benchmark iterations
-s                   show available input sets
-h                   this help text
Example:
To run 5 iterations with input set 10, execute "./md5 -i 10 -c 5"

## Inputs
Input sets are defined in the code (md5_bmark.c).
When added to this repository, the application offered the following input sets:
Index 0: 64 buffers of 512 bytes, totalling 0.0MB
Index 1: 64 buffers of 1024 bytes, totalling 0.1MB
Index 2: 64 buffers of 2048 bytes, totalling 0.1MB
Index 3: 64 buffers of 4096 bytes, totalling 0.3MB
Index 4: 128 buffers of 524288 bytes, totalling 67.1MB
Index 5: 128 buffers of 1048576 bytes, totalling 134.2MB
Index 6: 128 buffers of 2097152 bytes, totalling 268.4MB
Index 7: 128 buffers of 4194304 bytes, totalling 536.9MB
Index 8: 1000 buffers of 16384 bytes, totalling 16.4MB
Index 9: 2000 buffers of 16384 bytes, totalling 32.8MB
Index 10: 6400 buffers of 16384 bytes, totalling 104.9MB
Index 11: 1000 buffers of 32768 bytes, totalling 32.8MB
Index 12: 1000 buffers of 65536 bytes, totalling 65.5MB
Index 13: 1000 buffers of 131072 bytes, totalling 131.1MB
Index 14: 1000 buffers of 262144 bytes, totalling 262.1MB
Index 15: 1000 buffers of 524288 bytes, totalling 524.3MB
Index 16: 1000 buffers of 1048576 bytes, totalling 1048.6MB

## Additional information
###  Environment info
Vladimir (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7)
TaskSim tracing done with (MareNostrum): gcc (4.9.1) nanox (0.7a) mcc (1.99.7). For details ask Vladimir.

### Output with --summary
./md5 -i 16 -c 1
MSG: [?] ========== Nanos++ Initial Environment Summary ==========
MSG: [?] === PID:            10986
MSG: [?] === Num. threads:   16
MSG: [?] === Active CPUs:    [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, ]
MSG: [?] === Binding:        true
MSG: [?] === Prog. Model:    OmpSs
MSG: [?] === Plugin:         OpenCL PE Plugin
MSG: [?] ===  | Threads:     0
MSG: [?] =========================================================
MSG: [?] ============ Nanos++ Final Execution Summary ============
MSG: [?] === Application ended in 22 seconds
MSG: [?] === 1000 tasks have been executed
MSG: [?] =========================================================

