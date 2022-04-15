// ************************************************************************
//
// miniAMR: stencil computations with boundary exchange and AMR.
//
// Copyright (2014) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
// Questions? Contact Courtenay T. Vaughan (ctvaugh@sandia.gov)
//                    Richard F. Barrett (rfbarre@sandia.gov)
//
// ************************************************************************

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#if _POSIX_C_SOURCE < 199309L
#include <sys/time.h>
#endif

#include "block.h"
#include "proto.h"
#include "timer.h"

static inline double wtime(void)
{
#if _POSIX_C_SOURCE < 199309L
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (tv.tv_sec + 1e-6*tv.tv_usec);
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (ts.tv_sec + 1e-9*ts.tv_nsec);
#endif
}

double timer(void)
{
  // return((((double) clock())/((double) CLOCKS_PER_SEC)));
  return wtime();
}

void *ma_malloc(size_t size, char *file, int line)
{
   void *ptr;

   ptr = (void *) malloc(size);

   if (ptr == NULL) {
      printf("NULL pointer from malloc call in %s at %d\n", file, line);
      exit(-1);
   }

   counter_malloc++;
   size_malloc += (double) size;

   return(ptr);
}

void *ma_aligned_malloc(size_t size, char *file, int line)
{
  void *ptr = NULL;

  int res = posix_memalign(&ptr, sysconf(_SC_PAGESIZE), size);

  if (ptr == NULL) {
    printf("NULL pointer from posix_memalign call in %s at %d\n", file, line);
    exit(-1);
  }
  else if (res != 0) {
     printf("Unknown posix_memalign error (%d) in %s at %d\n", res, file, line);
     exit(-1);
  }

  counter_malloc++;
  size_malloc += (double) size;

  return(ptr);
}
