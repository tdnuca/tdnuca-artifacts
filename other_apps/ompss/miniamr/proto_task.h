#ifndef __PROTO_TASK_H
#define __PROTO_TASK_H

void init_task(size_t xdim, size_t ydim, size_t zdim,
               double *array,
               long int rand_seed);

void split_block_task(size_t xdim, size_t ydim, size_t zdim,
                      double *array_from,
                      double *array_to,
                      int i1, int j1, int k1);

void consolidate_block_task(size_t xdim, size_t ydim, size_t zdim,
                            double *array_from,
                            double *array_to,
                            int i1, int j1, int k1);

void stencil_task7(size_t xdim, size_t ydim, size_t zdim,
                   double *array);

void stencil_task27(size_t xdim, size_t ydim, size_t zdim,
                    double *array);

double check_sum_task(size_t xdim, size_t ydim, size_t zdim, double *array);

#endif
