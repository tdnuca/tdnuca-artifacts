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
#include "proto.h"
#include "proto_task.h"
#include <stdio.h>

// This routine does the stencil calculations.
void stencil_calc(int var, int number)
{
  int n, in;
  block *bp;

  if (stencil == 7) {
    for (in = 0; in < sorted_index[num_refine+1]; in++) {
      n = sorted_list[in].n;
      bp = &blocks[n];
      if (bp->number >= 0) {
        double * array = bp->array;
#pragma omp task label(stencil7)                                        \
  inout(([num_vars*box_size]array)[var*box_size;number*box_size])       \
  firstprivate(x_block_size, y_block_size, z_block_size,                \
               array, var, number)                                      \
  default(shared)
        {
          for (int i = var; i < var + number; ++i)
            stencil_task7(x_block_size+2,
                          y_block_size+2,
                          z_block_size+2,
                          &array[i*box_size]);
        }
      }
    }
  }
  else {
    for (in = 0; in < sorted_index[num_refine+1]; in++) {
      n = sorted_list[in].n;
      bp = &blocks[n];
      if (bp->number >= 0) {
        double * array = bp->array;
#pragma omp task label(stencil27)                                       \
  inout(([num_vars*box_size]array)[var*box_size;number*box_size])       \
  firstprivate(x_block_size, y_block_size, z_block_size,                \
               array, var, number)                                      \
  default(shared)
        {
          for (int i = var; i < var + number; ++i)
            stencil_task27(x_block_size+2,
                           y_block_size+2,
                           z_block_size+2,
                           &array[i*box_size]);
        }
      }
    }
  }
}

void stencil_calc_checksum(int var, int number, double *sum)
{
  int n, in;
  block *bp;

  if (stencil == 7) {
    for (in = 0; in < sorted_index[num_refine+1]; in++) {
      n = sorted_list[in].n;
      bp = &blocks[n];
      if (bp->number >= 0) {
        double * array = bp->array;
#pragma omp task label(stencil7_cs)                                     \
  inout(([num_vars*box_size]array)[var*box_size;number*box_size])       \
  concurrent(([num_vars]sum)[var;number])                               \
  firstprivate(x_block_size, y_block_size, z_block_size,                \
               array, var, number, sum)                                 \
  default(shared)
        for (int i = var; i < var + number; ++i) {
          stencil_task7(x_block_size+2,
                        y_block_size+2,
                        z_block_size+2,
                        &array[i*box_size]);
          
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
  }
  else {
    for (in = 0; in < sorted_index[num_refine+1]; in++) {
      n = sorted_list[in].n;
      bp = &blocks[n];
      if (bp->number >= 0) {
        double * array = bp->array;
#pragma omp task label(stencil27_cs)                                    \
  inout(([num_vars*box_size]array)[var*box_size;number*box_size])       \
  concurrent(([num_vars]sum)[var;number])                               \
  firstprivate(x_block_size, y_block_size, z_block_size,                \
               array, var, number)                                      \
  default(shared)
        for (int i = var; i < var + number; ++i) {
          stencil_task27(x_block_size+2,
                         y_block_size+2,
                         z_block_size+2,
                         &array[i*box_size]);
          
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
  }
}
