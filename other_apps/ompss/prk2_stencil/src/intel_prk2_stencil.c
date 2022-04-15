/*
Copyright (c) 2013, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
* Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

/*******************************************************************

NAME:    Stencil

PURPOSE: This program tests the efficiency with which a space-invariant,
         linear, symmetric filter (stencil) can be applied to a square
         grid or image.

USAGE:   The program takes as input the number of threads, the linear
         dimension of the grid, and the number of iterations on the grid

               <progname> <# threads> <iterations> <grid size>

         The output consists of diagnostics to make sure the
         algorithm worked, and of timing statistics.

FUNCTIONS CALLED:

         Other than OpenMP or standard C functions, the following
         functions are used in this program:

         wtime()
         bail_out()

HISTORY: - Written by Rob Van der Wijngaart, November 2006.
         - RvdW: Removed unrolling pragmas for clarity;
           added constant to array "in" at end of each iteration to force
           refreshing of neighbor data in parallel versions; August 2013

*******************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef ENABLE_MEMKIND
#include <memkind.h>
#endif

#ifdef ENABLE_PARSEC_HOOKS
#include "hooks.h"
#endif

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif
#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#endif

#define RESTRICT_KEYWORD

#ifdef RESTRICT_KEYWORD
#define RESTRICT restrict
#else
#define RESTRICT
#endif

#ifdef LONG_IS_64BITS
typedef unsigned long u64Int;
typedef long s64Int;
#define FSTR64 "%16ld"
#define FSTR64U "%16lu"
#else
typedef unsigned long long u64Int;
typedef long long s64Int;
#define FSTR64 "%16ll"
#define FSTR64U "%16llu"
#endif

#include <sys/time.h>
#define USEC_TO_SEC 1.0e-6 /* to convert microsecs to secs */
//extern double wtime(void);

#define DOUBLE

#define STAR

#ifndef RADIUS
  #define RADIUS 2
#endif

#define FULL_BLOCK_COPIES

#include "stencil_kernel.h"

#ifdef DOUBLE
  #define DTYPE   double
  #define EPSILON 1.e-8
  #define COEFX   1.0
  #define COEFY   1.0
  #define FSTR    "%lf"
#else
  #define DTYPE   float
  #define EPSILON 0.0001f
  #define COEFX   1.0f
  #define COEFY   1.0f
  #define FSTR    "%f"
#endif

void printBlock(int dim, int ts, int sj, int endj, int si, int endi, DTYPE (*m)[dim][dim], int phalo, int nofirstblock) {
   int startblock_i = sj + (nofirstblock ? ts : 0);
   int endblock_i   = endj - (nofirstblock ? ts : 0);
   printf("hola, %d %d (j %d, endj %d) (i %d, endi %d)\n", startblock_i , endblock_i, sj, endj, si, endi);
   for (int i = startblock_i; i < endblock_i; i += ts) {
      for ( int id = startblock_i; id < endblock_i; id+=ts ) printf("[ %d %p ]----------------------------------", i, &((*m)[i][id]) );
      printf("+\n");
      int start_i = phalo ? ( i == 0 ? RADIUS : i ) : i;
      int end_i   = phalo ? MIN(i+ts, dim-RADIUS) : i+ts;
      for (int ii = start_i; ii < end_i; ii += 1) {

         int startblock_j = si + (nofirstblock ? ts : 0);
         int endblock_j   = endi - (nofirstblock ? ts : 0);

         for (int j = startblock_j; j < endblock_j; j += ts) {
            int start_j = phalo ? ( j == 0 ? RADIUS : j ) : j;
            int end_j = phalo ? MIN(j+ts, dim-RADIUS) : j+ts;
            printf("%d | ", ii);
            for (int jj = start_j; jj < end_j; jj += 1) {
               printf("%.2f ", (*m)[ii][jj]);
            }
            printf("| ");
         }
         printf("\n");
      }
   }
}

void printMat(int dim, int ts, DTYPE (*m)[dim][dim], int phalo, int nofirstblock) {
   int startblock = nofirstblock ? ts : 0;
   int endblock   = nofirstblock ? dim-ts : dim;
   for (int i = startblock; i < endblock; i += ts) {
      for ( int id = startblock; id < endblock; id+=ts ) printf("[ %d %p ]----------------------------------", i, &((*m)[i][id]) );
      printf("+\n");
      int start_i = phalo ? ( i == 0 ? RADIUS : i ) : i;
      int end_i   = phalo ? MIN(i+ts, dim-RADIUS) : i+ts;
      for (int ii = start_i; ii < end_i; ii += 1) {


         for (int j = startblock; j < endblock; j += ts) {
            int start_j = phalo ? ( j == 0 ? RADIUS : j ) : j;
            int end_j = phalo ? MIN(j+ts, dim-RADIUS) : j+ts;
            for (int jj = start_j; jj < end_j; jj += 1) {
               printf("%.2f ", (*m)[ii][jj]);
            }
            printf("| ");
         }
         printf("\n");
      }
   }
}


double wtime() {
   double time_seconds;
   struct timeval time_data; /* seconds since 0 GMT */
   gettimeofday(&time_data,NULL);
   time_seconds = (double) time_data.tv_sec;
   time_seconds += (double) time_data.tv_usec * USEC_TO_SEC;
   return time_seconds;
}

/* define shorthand for indexing a multi-dimensional array                       */
#define IN(i,j)       in[i+(j)*(n)]
#define OUT(i,j)      out[i+(j)*(n)]
#define WEIGHT(ii,jj) (*weight)[ii+RADIUS][jj+RADIUS]

int main(int argc, char ** argv) {

  size_t    n;               /* linear grid dimension                            */
  long long int i, j, ii, jj, it, jt, iter;  /* dummies                          */
  DTYPE  norm, norm2,     /* L1 norm of solution                                 */
         reference_norm;
  DTYPE  f_active_points; /* interior of grid with respect to stencil            */
  DTYPE  f_active_points2;/* interior of grid with respect to stencil            */
  DTYPE  flops;           /* floating point ops per iteration                    */
  int    iterations;      /* number of times to run the algorithm                */
  double stencil_time,    /* timing parameters                                   */
         avgtime;
  int    stencil_size;    /* number of points in stencil                         */
  int    tile_size;       /* grid block factor                                   */
  size_t    total_length;    /* total required length to store grid values          */
  int    num_error=0;     /* flag that signals that requested and obtained
                             numbers of threads are the same                     */
  DTYPE  (* RESTRICT weight)[2*RADIUS+1][2*RADIUS+1]; /* weights of points in the stencil     */
  int radius = RADIUS;
  int use_hbm = 0;
  unsigned long total_alloc_size = 0;

  /*******************************************************************************
  ** process and test input parameters
  ********************************************************************************/

  if (argc != 5){
    printf("Usage: %s <# iterations> <array dimension> <tile size> [h|n] (h:use hbm n:do not use hbm)\n",
           *argv);
    return(EXIT_FAILURE);
  }

  iterations  = atoi(*++argv);
  if (iterations < 1){
    printf("ERROR: iterations must be >= 1 : %d \n",iterations);
    exit(EXIT_FAILURE);
  }

  n  = atoi(*++argv);

  if (n < 1){
    printf("ERROR: grid dimension must be positive: %lu\n", n);
    exit(EXIT_FAILURE);
  }

  if (RADIUS < 1) {
    printf("ERROR: Stencil radius %d should be positive\n", RADIUS);
    exit(EXIT_FAILURE);
  }

  if (2*RADIUS +1 > n) {
    printf("ERROR: Stencil radius %d exceeds grid size %lu\n", RADIUS, n);
    exit(EXIT_FAILURE);
  }

  DTYPE  (* RESTRICT in) [n][n]; /* input grid values                           */
  DTYPE  (* RESTRICT out)[n][n]; /* output grid values                          */

  /*  make sure the vector space can be represented                             */
  total_length = n*n*sizeof(DTYPE);
  if (total_length/n != n*sizeof(DTYPE)) {
    printf("ERROR: Space for %zu x %zu grid cannot be represented; [ %zu != %zu ]", n, n, total_length/n , n*sizeof(DTYPE));
    exit(EXIT_FAILURE);
  }

  tile_size = atoi(*++argv);
  if (tile_size < 1) {
     printf("ERROR: tile size must be positive : %d\n", tile_size);
     exit(EXIT_FAILURE);
  }
  char hbm_setting = (*++argv)[0];
  if (hbm_setting == 'h') {
     use_hbm = 1;
  } else if (hbm_setting == 'n') {
     use_hbm = 0;
  } else {
     printf("ERROR: hbm setting should be 'h' or 'n': %c\n", hbm_setting);
     exit(EXIT_FAILURE);
  }

  if ( use_hbm ) {
#ifdef ENABLE_MEMKIND
     printf("Allocating weight: %lu bytes\n",sizeof(*weight));
     total_alloc_size += sizeof(*weight);
     weight = (DTYPE (* RESTRICT)[2*RADIUS+1][2*RADIUS+1])memkind_malloc( MEMKIND_HBW_PREFERRED, sizeof(*weight) );

     printf("Allocating in: %.2lf GBs\n",(double)total_length/(1024.0*1024.0*1024.0));
     total_alloc_size += total_length;
     in  = (DTYPE (*)[n][n]) memkind_malloc(MEMKIND_HBW_PREFERRED, total_length);

     printf("Allocating out: %.2lf GBs\n",(double)total_length/(1024.0*1024.0*1024.0));
     total_alloc_size += total_length;
     out = (DTYPE (*)[n][n]) memkind_malloc(MEMKIND_HBW_PREFERRED, total_length);
#endif
  } else {
     printf("Allocating weight: %lu bytes\n",sizeof(*weight));
     total_alloc_size += sizeof(*weight);
     weight = (DTYPE (* RESTRICT)[2*RADIUS+1][2*RADIUS+1])malloc( sizeof(*weight) );

     printf("Allocating in: %.2lf GBs\n",(double)total_length/(1024.0*1024.0*1024.0));
     total_alloc_size += total_length;
     in  = (DTYPE (*)[n][n]) malloc(total_length);

     printf("Allocating out: %.2lf GBs\n",(double)total_length/(1024.0*1024.0*1024.0));
     total_alloc_size += total_length;
     out = (DTYPE (*)[n][n]) malloc(total_length);
  }
  if (!in || !out) {
    printf("ERROR: could not allocate space for input or output array\n");
    exit(EXIT_FAILURE);
  }

  /* fill the stencil weights to reflect a discrete divergence operator         */
  for (jj=-RADIUS; jj<=RADIUS; jj++) for (ii=-RADIUS; ii<=RADIUS; ii++)
    WEIGHT(ii,jj) = (DTYPE) 0.0;
#ifdef STAR
  stencil_size = 4*RADIUS+1;
  for (ii=1; ii<=RADIUS; ii++) {
    WEIGHT(0, ii) = WEIGHT( ii,0) =  (DTYPE) (1.0/(2.0*ii*RADIUS));
    WEIGHT(0,-ii) = WEIGHT(-ii,0) = -(DTYPE) (1.0/(2.0*ii*RADIUS));
  }
#else
  stencil_size = (2*RADIUS+1)*(2*RADIUS+1);
  for (jj=1; jj<=RADIUS; jj++) {
    for (ii=-jj+1; ii<jj; ii++) {
      WEIGHT(ii,jj)  =  (DTYPE) (1.0/(4.0*jj*(2.0*jj-1)*RADIUS));
      WEIGHT(ii,-jj) = -(DTYPE) (1.0/(4.0*jj*(2.0*jj-1)*RADIUS));
      WEIGHT(jj,ii)  =  (DTYPE) (1.0/(4.0*jj*(2.0*jj-1)*RADIUS));
      WEIGHT(-jj,ii) = -(DTYPE) (1.0/(4.0*jj*(2.0*jj-1)*RADIUS));
    }
    WEIGHT(jj,jj)    =  (DTYPE) (1.0/(4.0*jj*RADIUS));
    WEIGHT(-jj,-jj)  = -(DTYPE) (1.0/(4.0*jj*RADIUS));
  }
#endif

  norm = (DTYPE) 0.0;
  norm2 = (DTYPE) 0.0;
  f_active_points = (DTYPE) (n-2*RADIUS)*(DTYPE) (n-2*RADIUS);
  f_active_points2 = (DTYPE) (n-(2*tile_size))*(DTYPE) (n-(2*tile_size));

  printf("OmpSs stencil execution on 2D grid\n");
  printf("Grid size            = %lu\n", n);
  printf("Radius of stencil    = %d\n", RADIUS);
  if (tile_size <n-2*RADIUS)
     printf("Tile size            = %d\n", tile_size);
  else
     printf("Grid not tiled\n");
#ifdef STAR
  printf("Type of stencil      = star\n");
#else
  printf("Type of stencil      = compact\n");
#endif
#ifdef DOUBLE
  printf("Data type            = double precision\n");
#else
  printf("Data type            = single precision\n");
#endif
  printf("Number of iterations = %d\n", iterations);

#if 0
  /* intialize the input and output arrays                                     */
  for (j=0; j<n; j += tile_size) {
#pragma omp target device(smp) no_copy_deps copy_out( (*in)[j;tile_size][0;n] )
#pragma omp task label(top_init_in) out( (*in)[j][0] ) firstprivate(j, n, tile_size)
     {
        for (i=0; i<n; i += tile_size) {
//#pragma omp target device(smp) no_copy_deps copy_out( (*in)[j;tile_size][i;tile_size] )
//#pragma omp task label(init_in) out( (*in)[j][i] ) firstprivate(i, j, n, tile_size)
           for (jj=j; jj<j+tile_size; jj += 1) {
              for (ii=i; ii<i+tile_size; ii += 1) {
                 (*in)[jj][ii] = COEFX*((ii+RADIUS)-tile_size)+COEFY*((jj+RADIUS)-tile_size);
              }
           }
        }
     //#pragma omp taskwait noflush
     }
  }
  #pragma omp taskwait noflush

  for (j=tile_size; j<n-tile_size; j += tile_size) {
#pragma omp target device(smp) no_copy_deps copy_out( (*out)[j;tile_size][0;n] )
#pragma omp task label(top_init_out) out( (*out)[j][0] ) firstprivate(j, n, tile_size)
     {
        for (i=0; i<n; i += tile_size) {
//#pragma omp target device(smp) no_copy_deps copy_out( (*out)[j;tile_size][i;tile_size] )
//#pragma omp task label(init_out) out( (*out)[j][i] ) firstprivate(i, j, n, tile_size)
           {
              int first_j = (j == 0) ? RADIUS : j;
              int last_j = ((j+tile_size) >= n) ? n-RADIUS : j+tile_size;
              int first_i = (i == 0) ? RADIUS : i;
              int last_i = ((i+tile_size) >= n) ? n-RADIUS : i+tile_size;
              for (jj=first_j; jj<last_j; jj += 1) {
                 for (ii=first_i; ii<last_i; ii += 1) {
                    (*out)[jj][ii] = (DTYPE)0.0;
                 }
              }
           }
        }
//#pragma omp taskwait noflush
     }
  }
#else
  /* intialize the input and output arrays                                     */
  for (j=0; j<n; j += tile_size) {
#pragma omp target device(smp) no_copy_deps copy_out( (*in)[j;tile_size][0;n], (*out)[j;tile_size][0;n] )
#pragma omp task label(top_init_in) out( (*in)[j][0], (*out)[j][0] ) firstprivate(j, n, tile_size)
     {
        for (i=0; i<n; i += tile_size) {
           for (jj=j; jj<j+tile_size; jj += 1) {
              for (ii=i; ii<i+tile_size; ii += 1) {
                 (*in)[jj][ii] = COEFX*((ii+RADIUS)-tile_size)+COEFY*((jj+RADIUS)-tile_size);
              }
           }
        }
        for (i=0; i<n; i += tile_size) {
           {
              int first_j = (j == 0) ? RADIUS : j;
              int last_j = ((j+tile_size) >= n) ? n-RADIUS : j+tile_size;
              int first_i = (i == 0) ? RADIUS : i;
              int last_i = ((i+tile_size) >= n) ? n-RADIUS : i+tile_size;
              for (jj=first_j; jj<last_j; jj += 1) {
                 for (ii=first_i; ii<last_i; ii += 1) {
                    (*out)[jj][ii] = (DTYPE)0.0;
                 }
              }
           }
        }
     }
#pragma omp taskwait noflush
  }
#endif

  #pragma omp taskwait noflush

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_begin();
#endif

  for (iter = 0; iter<=iterations; iter++){

    /* start timer after a warmup iteration                                        */
    if (iter == 1) {
      #pragma omp taskwait noflush
      {
      //fprintf(stderr, "end of warm up.\n");
        stencil_time = wtime();
      }
    }

    /* Apply the stencil operator; only use tiling if the tile size is smaller
       than the iterior part of the grid                                       */
    if (tile_size < n-2*RADIUS) {
       for (j=tile_size; j<n-tile_size; j+=tile_size) {

          for (i=tile_size; i<n-tile_size; i+=tile_size) {
             DTYPE  (* RESTRICT in_up)    [n][n] = in;
             DTYPE  (* RESTRICT in_down)  [n][n] = in;
#ifdef FULL_BLOCK_COPIES
#pragma omp target device(smp) no_copy_deps \
             copy_inout( (*out)[j;tile_size][i;tile_size] ) \
             copy_in( *weight, (*in_up)    [j-tile_size;tile_size][i;tile_size], \
                               (*in)       [j          ;tile_size][i-tile_size;tile_size], \
                               (*in)       [j          ;tile_size][i          ;tile_size], \
                               (*in)       [j          ;tile_size][i+tile_size;tile_size], \
                               (*in_down)  [j+tile_size;tile_size][i;tile_size] )
#else
#pragma omp target device(smp) no_copy_deps \
             copy_inout( (*out)[j;tile_size][i;tile_size] ) \
             copy_in( *weight, (*in_up)    [j-radius;radius]   [i;tile_size], \
                               (*in)       [j       ;tile_size][i-radius;tile_size+2*radius], \
                               (*in_down)  [j+tile_size;radius][i;tile_size] )
#endif
#pragma omp task label (stencil) firstprivate(i, j, tile_size, radius, n, iter) \
             in( (*in)[j][i], (*in)[j][i+tile_size], (*in)[j][i-tile_size], (*in)[j-tile_size][i], (*in)[j+tile_size][i] )\
             inout( (*out)[j][i] )
                {
                    stencil_kernel_star_3in( n, in, in_up, in_down, out, j, i, tile_size, weight );
                }
            }
        }
    } else {
      //#pragma omp for
      for (j=RADIUS; j<n-RADIUS; j++) {
        for (i=RADIUS; i<n-RADIUS; i++) {
#ifdef STAR
          for (jj=-RADIUS; jj<=RADIUS; jj++)  (*out)[j][i] += WEIGHT(0,jj)*(*in)[j+jj][i];
          for (ii=-RADIUS; ii<0; ii++)        (*out)[j][i] += WEIGHT(ii,0)*(*in)[j][i+ii];
          for (ii=1; ii<=RADIUS; ii++)        (*out)[j][i] += WEIGHT(ii,0)*(*in)[j][i+ii];
#else
          /* would like to be able to unroll this loop, but compiler will ignore  */
          for (jj=-RADIUS; jj<=RADIUS; jj++)
          for (ii=-RADIUS; ii<=RADIUS; ii++)  (*out)[j][i] += WEIGHT(ii,jj)*(*in)[j+jj][i+ii];
#endif
        }
      }
    }

    /* add constant to solution to force refresh of neighbor data, if any       */
    for (j=0; j<n; j += tile_size) {
//#pragma omp target device(smp) copy_inout( (*in)[j;tile_size][0;n] ) no_copy_deps
////#pragma omp task label(update) inout( (*in)[j][i], (*in)[j+tile_size-1][i+tile_size-1], (*in)[j+tile_size-1][i], (*in)[j][i+tile_size-1]  ) firstprivate(n, tile_size, j, i)
//#pragma omp task label(top_update) inout( (*in)[j][0] ) firstprivate(n, tile_size, j)
       {
          //fprintf(stderr, "update w j=%d\n", j);
          for (i=0; i<n; i += tile_size) {
#pragma omp target device(smp) copy_inout( (*in)[j;tile_size][i;tile_size] ) no_copy_deps
#pragma omp task label(update) inout( (*in)[j][i] ) firstprivate(n, tile_size, j, i)
//#pragma omp target device(smp) no_copy_deps
//#pragma omp task label(update) priority(2)
             for (jj=j; jj<j+tile_size; jj++) {
                for (ii=i; ii<i+tile_size; ii++) {
                   (*in)[jj][ii] += 1.0;
                }
             }
          }
//#pragma omp taskwait noflush
       }
    }
  } /* end of iterations                                                        */

#ifdef VALIDATE
  #pragma omp taskwait
#else
  #pragma omp taskwait noflush
#endif
  {
    stencil_time = wtime() - stencil_time;
  }

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif

#ifdef VALIDATE
  // printMat(n, tile_size, in, 1, 1);
  // printf("#########################\n");
  // printMat(n, tile_size, out, 1, 1);
  /* compute L1 norm in parallel                                                */
  for (j=RADIUS; j<n-RADIUS; j++) for (i=RADIUS; i<n-RADIUS; i++) {
    norm += (DTYPE)ABS((*out)[j][i]);
  }

  for (j=tile_size; j<n-tile_size; j++) {
     for (i=tile_size; i<n-tile_size; i++) {
        norm2 += (DTYPE)ABS((*out)[j][i]);
     }
  }
  reference_norm = (DTYPE) (iterations+1) * (COEFX + COEFY);

  norm /= f_active_points;
  norm2 /= f_active_points2;
  printf("ref_norm=%f norm= %f, norm2=%f, f_active_points=%f, f_active_points2=%f time: %lf stencil_size=%d\n", reference_norm, norm, norm2, f_active_points, f_active_points2, stencil_time, stencil_size);

  /*******************************************************************************
  ** Analyze and output results.
  ********************************************************************************/

/* verify correctness                                                            */
  if (ABS(norm2-reference_norm) > EPSILON) {
    printf("ERROR: L1 norm = "FSTR", Reference L1 norm = "FSTR"\n",
           norm2, reference_norm);
    exit(EXIT_FAILURE);
  }
  else {
    printf("Solution validates\n");
#ifdef VERBOSE
    printf("Reference L1 norm = "FSTR", L1 norm = "FSTR"\n",
           reference_norm, norm2);
#endif
  }

#endif

  flops = (DTYPE) (2*stencil_size+1) * f_active_points2;
  avgtime = stencil_time/iterations;
  double gflops = 1.0E-09 * flops/avgtime;
  printf("Rate (MFlops/s): "FSTR"  Avg time (s): %lf, total_time: %lf\n",
         1.0E-06 * flops/avgtime, avgtime, stencil_time);
  // Print results
  printf("============= STENCIL RESULTS ============\n" );
  printf("  Execution time (sec): %f\n", stencil_time);
  printf("  Performance (GFLOPS): %f\n", gflops);
  printf("  Memory alocated (GB): %.2lf\n",(double)total_alloc_size/(1024.0*1024.0*1024.0));
  printf("==========================================\n" );

  exit(EXIT_SUCCESS);
}
