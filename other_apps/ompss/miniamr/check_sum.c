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

#include "block.h"
#include "timer.h"
#include "proto.h"
#include "proto_task.h"
#include <stdio.h>

// Generate check sum for a variable over all active blocks.
void check_sum(int var, int number, double *sum, int task)
{
  int n, in;
  double t1, t2;
  block *bp;

  t1 = timer();

  // memset((void*) &sum[var], 0, number*sizeof(double));
  for (in = 0; in < sorted_index[num_refine+1]; in++) {
    n = sorted_list[in].n;
    bp = &blocks[n];
    if (bp->number >= 0) {
      double * array = bp->array;
#pragma omp task label(check_sum_task)                            \
  concurrent(([num_vars]sum)[var;number])                         \
  inout(([num_vars*box_size]array)[var*box_size;number*box_size]) \
  firstprivate(x_block_size, y_block_size, z_block_size,          \
               number, var, sum, array, n, in)                    \
  if(task)
      for (int i = var; i < var + number; ++i) {
         double output =
           check_sum_task(x_block_size+2,
                          y_block_size+2,
                          z_block_size+2,
                          &array[i*box_size]);

#pragma omp atomic
         sum[i] += output;
      }
    }
  }

  t2 = timer();
  timer_cs_calc += t2 - t1;
  total_red++;
}
