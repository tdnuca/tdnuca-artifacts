# Integral histogram (crossweave scan)

Author: Pieter Bellens (BAR)

Modified by Isaac Sánchez Barrera <isaac.sanchez@bsc.es>

See [1, 2].


## Compilation instructions

```
autoreconf -fi
./configure --prefix=out
make
make install
```

This will install the binaries in `$PWD/out/bin`.

There are some interesting flags for the `configure` step:
- `--disable-ompss`: To disable OmpSs calls and generate a sequential app.
- `--disable-sequence`: To disable *sequence mode* in the ompss version. The
  sequence mode has no taskwaits between iterations and allocates twice the
  amount of memory for the halos and histograms (similar to a  double buffer).
  Only for OmpSs.
- `--disable-socket`: To disable the manual socket-aware calls. Only for OmpSs.

And there are some useful variables, too:
- `CC`: For the C compiler. If in OmpSs mode, the default compilers in order are
  `smpcc`, `mcc`, `imcc` and `xlmcc`. In sequential mode, the default compilers
  are `gcc`, `icc`, `xlc` and `cc`. Remember to use `--disable-ompss` if you
  want to use a non-OmpSs compiler, the configure step will fail otherwise.
- `CFLAGS`: The C compiler flags. By default, `-Wall`. For OmpSs, `--ompss` is
  added automatically.
- `CFLAGS_PERF`: Extra flags for the performance binary. By default, `-O3`.
- `CFLAGS_INSTR`: Extra flags for the instrumentation binary. By default, `-O3`;
  for OmpSs, `--instrument` is added automatically.
- `CFLAGS_DEBUG`: Extra flags for the instrumentation binary. By default, `-g`;
  for OmpSs, `--debug -k` are added automatically.
- `LIBS`: Extra linking libraries. By default, `-lm`. If you change it, do not
  forget to include `-lm` or equivalent.
- Similarly, there are `LDFLAGS`, and the `PERF`, `INSTR`, `DEBUG` variants.


## Running instructions and inputs

From `orig_README` (updated here).

The algorithm initializes each pixel with an empty histogram. The
number of bins and the size of a bin are input parameters. The
application computes the cumulative histogram at a pixel (i,j) in two
passes over the blocks of the input image. The horizontal pass
processes the image blocks and propagates the histograms from (i,j) to
(i,j+1) The vertical pass propagates the histograms from (i,j) to
(i+1,j).

`inthist_crossweave imw imh bw bh bins binstride ims [output]`

Where
- `imw`: image width in pixels
- `imh`:  image height in pixels
- `bw`:  width of a block in pixels
- `bh`:  height of a block in pixels
- `bins`: number of bins in a histogram
- `binstride`: bin width (floating point)
- `ims`: number of iterations (images)
- `output`: filename for an output of the timing stats. Optional.

To check the results for correctness define `INTHIST_CHECK` to a value greater
than 0 (either in `CFLAGS` or by editing `inthist_extra.h`).

The instrumentation and debug binaries are `inthist_crossweave_i` and 
`inthist_crossweave_d` respectively.


## Additional information

- Tested with the following (Isaac):
  * Mercurium 2.0.0 (abc233dd, 14th Oct 2016) and Nanox++ 0.10
    (3550efdd, 16th Jun 2016). Backend GCC 5.1.0 (MareNostrum)
- The code is NUMA-aware when using Nanos++ with
  `--schedule=socket --no-socket-auto-detect` (if not disabled during the
  configure step).
- All the matrix blocks are aligned to the page size, obtained by calling
  `sysconf(_SC_PAGESIZE)`. You can change it to value `size` by adding
  `-D_PAGE_SIZE=size` to the `CFLAGS`.
- Removed OpenCV and CUDA support.
 

### References

[1] Porikli, Fatih. "Integral histogram: a fast way to extract histograms in
	Cartesian spaces," _CVPR_, Sep. 2005. DOI: `<10.1109/CVPR.2005.188>`

[2] <https://pm.bsc.es/projects/bar/wiki/inthist_crossweave>


### Output using `--summary`


```
MSG: [?] ========== Nanos++ Initial Environment Summary ==========
MSG: [?] === PID:                 17430
MSG: [?] === Num. worker threads: 16
MSG: [?] === System CPUs:         active[ 0-15 ] - inactive[  ]
MSG: [?] === Binding:             true
MSG: [?] === Prog. Model:         OmpSs
MSG: [?] === Priorities:          Not needed / disabled
MSG: [?] === Plugin:              SMP PE Plugin
MSG: [?] ===  | PEs:              16
MSG: [?] ===  | Worker Threads:   16
MSG: [?] =========================================================
image 8x8 blocks 32x32 binsize 4 binstride 1.000000 imcnt 4
inthist wall clock            : 5338 us
inthist wall clock per image  : 1334.50 us
inthist frate                 : 749.34 fps
inthist blockrate             : 47958.04 blps
ERROR = 0
MSG: [?] ============ Nanos++ Final Execution Summary ============
MSG: [?] === Application ended in 0 seconds
MSG: [?] === 992 tasks have been executed
MSG: [?] =========================================================
```
