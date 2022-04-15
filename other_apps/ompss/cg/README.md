# A simple implementation of CG
parallelized with [OmpSs](https://pm.bsc.es/ompss) by Luc Jaulmes <luc.jaulmes@bsc.es>

### Ompss features

This code uses the following OmpSs features:

* Annotations of all dynamic dependencies.
  * Caveat: the matrix never changes, hence it is not really a dependency. However, it is the biggest part of the data
and is used (only) in tasks `compute_Ap` and `recompute_gradient_mvm`. To add dependencies on the matrix, you can
add for each block (s:e) the following: `A->r[s:e], A->c[A->r[s]:A->r[e]], A->v[A->r[s]:A->r[e]]`
* Multi-dependencies (see above)
  * This allows to avoid sentinels and keep simple dependency trackers
* Variable input sizes, see [below](#running). Two efficient paramters are to impact input size and
running time are the choice of the matrix and the maximum number of iterations.
* Overlapping of task creation with task execution (i.e. use of `taskwait on(...)`).
  * It is recommended (for the 16.06 release of OmpSs) to use `--disable-immediate-successor` in `NX_ARGS`.
  * A useful trick can be to use N-1 blocks on N threads to improve parallelism, dedicating one thread fully to
runtime tasks such as task creation.
* Task reductions with overlap (passing the `--enable-reductions` option to configure, this is not yet supported in the 16.06 release)
* Use of `concurrent(...)` for reductions while task reductions are unavailable
* Hybrid MPI+OmpSs parallelism (passing the `--with-mpi` option to configure)
  * A parallel efficiency (strong scaling) of 80.17% has been achieved with  matrix `-synth Poisson3D 27 512` on
128 sockets (1024 cores), with baseline 8 sockets (64 cores, smallest amount to fit the matrix in memory),
running on Intel Xeon E5-2670 [[Jaulmes et al., SC15]](https://dl.acm.org/citation.cfm?id=2807599).

It does not feature:

* Task nesting
* Out of the box OpenMP4.0 compatibility

### Compiling

To compile, you need         |  configure variable to adjust
------|------
a C compiler that understands c11 | `CC`
the mercurium compiler for OmpSs  | `MCC`
nanos                             | `NANOS_HOME`
extrae (optional)                 | `EXTRAE_HOME`
mpi (optional)                    | `MPI_HOME`

By default, 3 targets are available:
* `cg`, parallel version that is compiled with the performance library of nanos
* `cg_instr`, parallel version that is compiled with the instrumentation library of nanos
* `cg_seq`, sequential version compiled ignoring pragmas (entirely compiled with `$(CC)`)

Compile-time options that can be fed to the configure script are
* *verbosity level*: silent or varying from 1 (a little talkative) to 4 (way too much), see
`./configure --enable-verbosity=[level]`
* *optimization options*: The easiest way is to set these in `CFLAGS`.
* *extrae events*: Performance measures (solving start/end, convergence value and time per iteration)
can be written to a trace file instead of printing to stdout. Can be set with `--enable-extrae`
* *MPI*: Builds hybrid MPI+OmpSs, disable with `--without-mpi`, enable with `--with-mpi[=path]`,
optionally setting the path where MPI is installed this way.
* *task reductions*: To enable task reduction pragmas (instead of atomic concurrent accesses)
use `--enable-reductions` and to disable `--disable-reductions`.
Using reductions conflicts with overlapping task creation and computation, whic CG does.
This bug is fixed in nanox >= Jul 05, 2016 (commit 7680b4c), which is not included in the 16.06 release of OmpSs.

### Running

You have to supply a (sparse) **symmetric positive definite (spd) matrix** in one of two formats:
* a file containing the matrix in Matrix Market format (*not for MPI+OmpSs*). For inputs, see for example:<br />
  https://www.cise.ufl.edu/research/sparse/matrices/<br />
  http://math.nist.gov/MatrixMarket/
* a set of parameters for the program to generate one dynamically.<br />
  Only 3D Poisson's equation using finite differences is supported so far, with a 7, 19 or 27 points stencil.
  You can specify this with `-synth Poisson3D p n`, with p one of 7, 19, 27 and n^3 the total matrix size.

All other parameters are optional and described below.

```
Usage: ./cg [options] <matrix> [, ...]
Possible options are:
  -th    threads    Manually define number of threads.
  -nb    blocks     Defines the number of blocks in which to divide operations; blocks' actual
                    sizes will depend on the matrix' size. Usually equal to number of threads.
  -r     runs       number of times to run a matrix solving.
  -cv    thres      Run until the error verifies ||b-Ax|| < thres * ||b|| (default 1e-10).
  -maxit N          Run no more than N iterations (default no limit).
  -seed  s          Initialize seed of each run with s. If 0 use different (random) seeds.
All options apply to every following input file. You may re-specify them for each file.
```

e.g. the following command will run CG twice with the four blocks on four threads, and twice
with 8 blocks:

```
$ ./cg -r 2 -th 4 -nb 4 /path/to/file.mtx -nb 8 -synth Poisson3D 7 5
```

To run `cg` or `cg_seq` with extrae events enabled, set `EXTRAE_CONFIG_FILE` to [extrae\_conf.xml](../master/scripts/extrae_conf.xml)
This is a minimal configuration file for Extrae. With `cg_instr`, you may use this file or another
if you plan to use nanos' built-in Extrae instrumentation.

### Viewing

In the subdirectory [scripts/](../master/scripts/) you will find paraver configuration files that allow to view
extrae events that are specific to this application.

The [trace\_to\_plot.py](../master/scripts/trace_to_plot.py) script extracts these informations from configuration files into a
format that can be fed to any viewing tool, to plot convergence over time or iterations.

