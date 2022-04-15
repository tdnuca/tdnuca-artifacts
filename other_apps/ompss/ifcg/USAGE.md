# Iteration-Fusing Conjugate Gradient (IFCG)
Author: Sicong Zhuang <sicong.zhuang@bsc.es>

## Issues
Using the #pragma omp register clause introduced since the OmpSs 14.04.
Any lower versions are not supported.

## Compilation instructions
`make`

## Running instructions
USAGE: ./dcg [bm] [it] [precision] [correction] [iter_fuse] [rep] [orth_fac] [HB_MAT] [FULL?] [CG_VER] [cglog] [B]

bm : block size
it : maximun iterations
precision : desired precision
correction : accuracy improvement (suggested to be the same as iter_fuse if iter_fuse>1)
iter_fuse : number of fusing iterations
rep  : repetitions
orth_fac : 1E5 (no real effect in the current version)
HB_MAT : input matrix
FULL? : if the input matrix is full or truncated by the diagonal.
CG_VER : version of CGs
		0 : Preconditioned CG
		1 : Chronopoulos CG [1]
		2 : Pipelined CG [2]
		3 : Gropp CG [3]
		4 : IFCG
		5 : IFCG2
cglog : log file
B : right hand side file (if not specified, then will be randomly generated)

## Inputs
Requires a real-valued symmetric positive definite sparse matrix as the input in the form of a Harwell-Boeing or its variant Rutherford-Boeing format. 

A collection of such matrices is available at The SuiteSparse Matrix Collection
(http://www.cise.ufl.edu/research/sparse/matrices).

The current matrix reading routine imposes a strict format on the input matrix which the matrices from the aforementioned collection might not meet.
src/script/hb_parser.sh is an auxilary script to convert the meta-data.

An example of such matrix can be found at src/input/LFAT5.rb

## Additional information
Working under nanox 0.12a, mcxx 2.0.0, gcc 6.3.1

./dcg 2066 20 1E-8 10 10 5 1E5 qa8fm/qa8fm.rb 0 4 1
MSG: [?] ========== Nanos++ Initial Environment Summary ==========
MSG: [?] === PID:                 4575
MSG: [?] === Num. worker threads: 4
MSG: [?] === System CPUs:         active[ 0-3 ] - inactive[  ]
MSG: [?] === Binding:             true
MSG: [?] === Prog. Model:         OmpSs
MSG: [?] === Priorities:          Needed / enabled
MSG: [?] === Plugin:              SMP PE Plugin
MSG: [?] ===  | PEs:              4
MSG: [?] ===  | Worker Threads:   4
MSG: [?] =========================================================
n 66127 bm 2066 cgit 20 prec 1.000000E-08 correction 10 iter_fuse 10 rep 5 orth 1.000000E+05 CG_VER 4 A qa8fm/qa8fm.rb B (null) loglevel 1
gen rhs
warning: writing obj RHS
ALG4 IFCG
ALG4 IFCG
ALG4 IFCG
ALG4 IFCG
ALG4 IFCG
MSG: [?] ============ Nanos++ Final Execution Summary ============
MSG: [?] === Application ended in 5 seconds
MSG: [?] === 21650 tasks have been executed
MSG: [?] =========================================================

## References:
[1]A. T. Chronopoulos. 1991. s-Step Iterative Methods for (Non)Symmetric
(In)Denite Linear Systems. SIAM J. Numer. Anal. 28(6) (1991), 1776–1789.

[2]P.Ghysels and W.Vanroose. 2014. Hiding global synchronization latency in
the preconditioned Conjugate Gradient algorithm. Parallel Comput. 40 (2014),
224–238.

[3]W. Gropp. 2010. Update on libraries for blue waters. hp://jointlab-pc.ncsa.
illinois.edu/events/workshop3/pdf/presentations/Gropp-Update-on-Libraries.
pdf. (2010).
