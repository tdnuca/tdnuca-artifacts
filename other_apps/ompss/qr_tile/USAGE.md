# Tiled QR factorisation (C++)

Author: Isaac Sánchez Barrera <isaac.sanchez@bsc.es>

See Buttari et al. [1, 2].


## Compilation instructions

Edit the Makefile and set `LIBS` according to the path to OpenBLAS, MKL, or
whatever implementation of LAPACK you are using. Then just call `make` and
multiple binaries will be built.

**Important**: Make sure to link with the _sequential_ (read: not parallel)
version of the LAPACK implementation you use.


## Running instructions and inputs

There are two ways of executing the program, using a random matrix or reading
it from an input file.

For a random matrix, use `qr <block size> <total rows> <total columns>`. You
can add two extra parameters to save the Q and R matrices, although they are
useless in this case.

For an input matrix, use `qr <block size> <input file>`. The input file must
be in plain text format and start with two integers (the M rows and N columns)
followed by as many floating point numbers as MxN. For example:

```
3 4
1.0 1.0 1.0  1.0
1.0 2.0 4.0  8.0
1.0 3.0 9.0 27.0

```

Not too tested, but "good" block sizes for MareNostrum seem to have between
256 and 512 elements per dimension.


## Additional information

- Tested with the following (Isaac):
  * Mercurium 2.0.0 (abc233dd, 14th Oct 2016) and Nanox++ 0.10
    (3550efdd, 16th Jun 2016). Backend GCC 5.1.0 (MareNostrum)
  * LAPACK implementations:
    - OpenBLAS 0.2.15 (compiled with USE_THREADS=0)
    - Intel MKL (using `module load mkl` in MN)
- The code is NUMA-aware when using Nanos++ with
  `--schedule=socket --no-socket-auto-detect`.
- All the matrix blocks are aligned to the page size, obtained by calling
  `sysconf(_SC_PAGESIZE)`. You can change it to value `size` by adding
  `-D_PAGE_SIZE=size` to the `CXXFLAGS` in the Makefile.
- Although you need to set the block size for Nanos tasks, the program makes a
  call to the LAPACK function `ilaenv` to set the internal block size for
  LAPACK to a sensible value (which is implementation and architecture
  dependent). The NANOS tasks still have the granularity of the command line
  block size, but LAPACK internally subdivides the block into smaller blocks
  to achieve better performance. Should you wish to use a custom block size
  for LAPACK, change the value of `QR::LAPACK_bs`.
- The computation tasks (`QR::dgeqrt`, `QR::dlarfb`, `QR::dtpqrt` and
  `QR::dtpmqrt`) allocate and deallocate a `work` array used by the LAPACK
  procedures. They must be able to hold at least `BS*LAPACK_bs` double
  precision floating point values.
- If the numbers of rows or columns are not a multiple of the task block size,
  the matrix is grown until they are and the extra elements are set to zero.
  This does not change the correctness of the algorithm.
 

### References

[1] Buttari, Alfredo; Langou, Julien; Kurzak, Jakub, and Dongarra, Jack.
    "Parallel tiled QR factorization for multicore architetures,"
    _Concurrency Computat.: Pract. Exper._, vol. 20, no. 13, pp. 1573-1590,
    Sep. 2008. DOI: `<10.1002/cpe.1301>`

[2] Buttari, Alfredo; Langou, Julien; Kurzak, Jakub, and Dongarra, Jack.
    "A class of parallel tiled linear algebra algorithms for multicore
    architectures," _Parallel Computing_, vol. 35, no. 1, pp. 38-53,
    Jan. 2009. DOI: `<10.1016/j.parco.2008.10.002>`


### Output using `--summary`

For `qr 256 10240 10240` in MareNostrum using OpenBLAS:

```
MSG: [?] ========== Nanos++ Initial Environment Summary ==========
MSG: [?] === PID:                 15236
MSG: [?] === Num. worker threads: 16
MSG: [?] === System CPUs:         active[ 0-15 ] - inactive[  ]
MSG: [?] === Binding:             true
MSG: [?] === Prog. Model:         OmpSs
MSG: [?] === Priorities:          Not needed / disabled
MSG: [?] === Plugin:              SMP PE Plugin
MSG: [?] ===  | PEs:              16
MSG: [?] ===  | Worker Threads:   16
MSG: [?] =========================================================
Time: 5622635
LAPACK block size: 32
       1431.656 GFlop => 254623.64 MFlop/s (standard QR)
       1454.025 GFlop => 258602.13 MFlop/s (true calculations)
MSG: [?] ============ Nanos++ Final Execution Summary ============
MSG: [?] === Application ended in 6 seconds
MSG: [?] === 3670 tasks have been executed
MSG: [?] =========================================================
```
