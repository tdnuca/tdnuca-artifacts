#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// things that should be defined globally : constants, functions, etc.
#ifndef RECOMPUTE_GRADIENT_FREQ
#define RECOMPUTE_GRADIENT_FREQ 50
#endif

#ifdef UNUSED
#elif defined(__GNUC__)
	#define UNUSED __attribute__((unused))
#elif defined(__LCLINT__)
	#define UNUSED /*@unused@*/
#elif defined(__cplusplus)
	#define UNUSED
#endif

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>

// a few global vars for parameters that everyone needs to know
#ifndef USE_MPI
extern int nb_blocks;
extern int max_it;
extern int *block_bounds;

static inline void set_block_end(const int b, const int pos)
{
	block_bounds[b+1] = pos;
}
#else
extern int max_it;
// mpi_zone{start,size} gives repartition of matrix rows among MPI ranks
// block_bounds gives repartition of matrix rows among blocks (~ threads) within MPI rank
// glob_blocks gives global access indexes for local blocks,
// plus two globally-indexed blocks for non-local data (one before and one after local data)
extern int nb_blocks, mpi_rank, mpi_size, nb_glob_blocks;
extern int *block_bounds, *mpi_zonestart, *mpi_zonesize, *glob_block_bounds;
#endif

static inline int get_block_start(const int b)
{
	return block_bounds[b];
}

static inline int get_block_end(const int b)
{
	return block_bounds[b+1];
}

#ifdef USE_MPI
static inline int world_block(const int b)
{
	return nb_blocks * mpi_rank + b;
}

static inline double* localrank_ptr(double *local_ptr)
{
       return local_ptr + mpi_zonestart[mpi_rank];
}

static inline int glob_block_start(const int b)
{
	return glob_block_bounds[b];
}

static inline int glob_block_end(const int b)
{
	return glob_block_bounds[b+1];
}
#else
// shorthands so we can use the same function with and without MPI
static inline int world_block(const int b) {return b;}
static inline double *localrank_ptr(double *p) {return p;}
#endif

// these are defined here so that they are compiled with mcc
#ifndef _OMPSS
static inline int get_num_threads() { return 1; }
static inline int get_thread_num()  { return 0; }
static inline void set_num_threads(int t UNUSED)
{
	fprintf(stderr, "In sequential CG, num_threads can't be set, it is fixed to 1.\n");
	exit(2);
}
#else
#include <nanox/nanos_omp.h>
#define get_num_threads  nanos_omp_get_num_threads
#define get_thread_num   nanos_omp_get_thread_num
#define set_num_threads  nanos_omp_set_num_threads
#endif

#endif // GLOBAL_H_INCLUDED

