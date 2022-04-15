#include <stdlib.h>
#include <stdio.h>
#include "block.h"
#include "proto.h"

static inline void comm_eq(double *from,
                           double *to,
                           int itot, int iinc,
                           int jtot, int jinc)
{
  for (int i = 0; i < itot*iinc; i += iinc)
    for (int j = 0; j < jtot*jinc; j += jinc)
      to[i + j] = from[i + j];
}

static inline void comm_finetocoarse(double *from,
                                     double *to,
                                     int itot, int iinc,
                                     int jtot, int jinc)
{
  for (int i = 0; i < itot*iinc; i += iinc)
    for (int j = 0; j < jtot*jinc; j += jinc)
      to[i + j] =
        (from[i*2 + j*2] +
         from[i*2 + iinc + j*2 + jinc]) +
        (from[i*2 + j*2 + jinc] +
         from[i*2 + iinc + j*2]);
}

static inline void comm_coarsetofine(double *from,
                                     double *to,
                                     int itot, int iinc,
                                     int jtot, int jinc)
{
  for (int i = 0; i < itot*iinc; i += iinc)
    for (int j = 0; j < jtot*jinc; j += jinc)
      to[i*2 + j*2] =
        to[i*2 + j*2 + jinc] =
        to[i*2 + iinc + j*2] =
        to[i*2 + iinc + j*2 + jinc] =
        from[i + j]/4.0;
}

void comm_alt(double * east[4], int east_f,
              double * west[4], int west_f,
              double * north[4], int north_f,
              double * south[4], int south_f,
              double * up[4], int up_f,
              double * down[4], int down_f,
              double * b_array,
              int start, int number)
{
  int east_c = east_f == 1 ? 4 : 1;
  int west_c = west_f == 1 ? 4 : 1;
  int north_c = north_f == 1 ? 4 : 1;
  int south_c = south_f == 1 ? 4 : 1;
  int up_c = up_f == 1 ? 4 : 1;
  int down_c = down_f == 1 ? 4 : 1;
  
#pragma omp task label(comm_alt)                                        \
  inout(([num_vars*box_size]b_array)[start*box_size;number*box_size])   \
  in({([num_vars*box_size](east[i]))[start*box_size;number*box_size] , i=0;east_c}, \
     {([num_vars*box_size](west[i]))[start*box_size;number*box_size] , i=0;west_c}, \
     {([num_vars*box_size](north[i]))[start*box_size;number*box_size] , i=0;north_c}, \
     {([num_vars*box_size](south[i]))[start*box_size;number*box_size] , i=0;south_c}, \
     {([num_vars*box_size](up[i]))[start*box_size;number*box_size] , i=0;up_c}, \
     {([num_vars*box_size](down[i]))[start*box_size;number*box_size] , i=0;down_c}) \
  firstprivate(start, number,                                           \
               east_f, west_f, north_f, south_f, up_f, down_f,          \
               east_c, west_c, north_c, south_c, up_c, down_c)
  {
    typedef double (*box_t)[y_block_size+2][z_block_size+2];
    // EAST -> CENTER
    if (east_f == 2) {
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(east[0])[var*box_size])[x_block_size][1][1],
                &((box_t) &b_array[var*box_size])[x_block_size+1][1][1],
                y_block_size, z_block_size + 2,
                z_block_size, 1);
      }
    else if (east_f == 0) {
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(east[0])[var*box_size])[1][1][1],
                &((box_t) &b_array[var*box_size])[x_block_size+1][1][1],
                y_block_size, z_block_size + 2,
                z_block_size, 1);
    }
    else if (east_f < 0) {
      east_f = -east_f - 1;
      int ystart = 1 + (east_f%2)*y_block_half;
      int zstart = 1 + (east_f/2)*z_block_half;
      for (int var = start; var < start + number; ++var)
        comm_coarsetofine(&((box_t) &(east[0])[var*box_size])[1][ystart][zstart],
                          &((box_t) &b_array[var*box_size])[x_block_size+1][1][1],
                          y_block_half, z_block_size + 2,
                          z_block_half, 1);
    }
    else {      
      for (int var = 0; var < start + number; ++var)
        for (int i = 0; i < 4; ++i) {
          int ystart = 1 + (i%2)*y_block_half;
          int zstart = 1 + (i/2)*z_block_half;
          comm_finetocoarse(&((box_t) &(east[i])[var*box_size])[1][1][1],
                            &((box_t) &b_array[var*box_size])[x_block_size+1][ystart][zstart],
                            y_block_half, z_block_size + 2,
                            z_block_half, 1);
        }
    }

    fflush(stdout);
    
    // WEST -> CENTER
    if (west_f == 2) {
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(west[0])[var*box_size])[1][1][1],
                &((box_t) &b_array[var*box_size])[0][1][1],
                y_block_size, z_block_size + 2,
                z_block_size, 1);
    }
    else if (west_f == 0) {
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(west[0])[var*box_size])[x_block_size][1][1],
                &((box_t) &b_array[var*box_size])[0][1][1],
                y_block_size, z_block_size + 2,
                z_block_size, 1);
    }
    else if (west_f < 0) {
      west_f = -west_f - 1;
      int ystart = 1 + (west_f%2)*y_block_half;
      int zstart = 1 + (west_f/2)*z_block_half;
      for (int var = start; var < start + number; ++var)
        comm_coarsetofine(&((box_t) &(west[0])[var*box_size])[x_block_size][ystart][zstart],
                          &((box_t) &b_array[var*box_size])[0][1][1],
                          y_block_half, z_block_size + 2,
                          z_block_half, 1);
    }
    else {
      for (int var = 0; var < start + number; ++var)
        for (int i = 0; i < 4; ++i) {
          int ystart = 1 + (i%2)*y_block_half;
          int zstart = 1 + (i/2)*z_block_half;
          comm_finetocoarse(&((box_t) &(west[i])[var*box_size])[x_block_size][1][1],
                            &((box_t) &b_array[var*box_size])[0][ystart][zstart],
                            y_block_half, z_block_size + 2,
                            z_block_half, 1);
        }
    }

    fflush(stdout);

    // NORTH -> CENTER
    if (north_f == 2)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(north[0])[var*box_size])[1][y_block_size][1],
                &((box_t) &b_array[var*box_size])[1][y_block_size+1][1],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                z_block_size, 1);
    else if (north_f == 0)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(north[0])[var*box_size])[1][1][1],
                &((box_t) &b_array[var*box_size])[1][y_block_size+1][1],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                z_block_size, 1);
    else if (north_f < 0) {
      north_f = -north_f - 1;
      int xstart = 1 + (north_f%2)*x_block_half;
      int zstart = 1 + (north_f/2)*z_block_half;
      for (int var = start; var < start + number; ++var)
        comm_coarsetofine(&((box_t) &(north[0])[var*box_size])[xstart][1][zstart],
                          &((box_t) &b_array[var*box_size])[1][y_block_size+1][1],
                          x_block_half, (z_block_size + 2)*(y_block_size + 2),
                          z_block_half, 1);
    }
    else
      for (int var = 0; var < start + number; ++var)
        for (int i = 0; i < 4; ++i) {
          int xstart = 1 + (i%2)*x_block_half;
          int zstart = 1 + (i/2)*z_block_half;
          comm_finetocoarse(&((box_t) &(north[i])[var*box_size])[1][1][1],
                            &((box_t) &b_array[var*box_size])[xstart][y_block_size+1][zstart],
                            x_block_half, (z_block_size + 2)*(y_block_size + 2),
                            z_block_half, 1);
        }

    // SOUTH -> CENTER
    if (south_f == 2)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(south[0])[var*box_size])[1][1][1],
                &((box_t) &b_array[var*box_size])[1][0][1],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                z_block_size, 1);
    else if (south_f == 0)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(south[0])[var*box_size])[1][y_block_size][1],
                &((box_t) &b_array[var*box_size])[1][0][1],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                z_block_size, 1);
    else if (south_f < 0) {
      south_f = -south_f - 1;
      int xstart = 1 + (south_f%2)*x_block_half;
      int zstart = 1 + (south_f/2)*z_block_half;
      for (int var = start; var < start + number; ++var)
        comm_coarsetofine(&((box_t) &(south[0])[var*box_size])[xstart][y_block_size][zstart],
                          &((box_t) &b_array[var*box_size])[1][0][1],
                          x_block_half, (z_block_size + 2)*(y_block_size + 2),
                          z_block_half, 1);
    }
    else
      for (int var = 0; var < start + number; ++var)
        for (int i = 0; i < 4; ++i) {
          int xstart = 1 + (i%2)*x_block_half;
          int zstart = 1 + (i/2)*z_block_half;
          comm_finetocoarse(&((box_t) &(south[i])[var*box_size])[1][y_block_size][1],
                            &((box_t) &b_array[var*box_size])[xstart][0][zstart],
                            x_block_half, (z_block_size + 2)*(y_block_size + 2),
                            z_block_half, 1);
        }
    
    // UP -> CENTER
    if (up_f == 2)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(up[0])[var*box_size])[1][1][z_block_size],
                &((box_t) &b_array[var*box_size])[1][1][z_block_size+1],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                y_block_size, z_block_size + 2);
    else if (up_f == 0)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(up[0])[var*box_size])[1][1][1],
                &((box_t) &b_array[var*box_size])[1][1][z_block_size+1],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                y_block_size, z_block_size + 2);
    else if (up_f < 0) {
      up_f = -up_f - 1;
      int xstart = 1 + (up_f%2)*x_block_half;
      int ystart = 1 + (up_f/2)*y_block_half;
      for (int var = start; var < start + number; ++var)
        comm_coarsetofine(&((box_t) &(up[0])[var*box_size])[xstart][ystart][1],
                          &((box_t) &b_array[var*box_size])[1][1][z_block_size+1],
                          x_block_half, (z_block_size + 2)*(y_block_size + 2),
                          y_block_half, z_block_size + 2);
    }
    else
      for (int var = 0; var < start + number; ++var)
        for (int i = 0; i < 4; ++i) {
          int xstart = 1 + (i%2)*x_block_half;
          int ystart = 1 + (i/2)*y_block_half;
          comm_finetocoarse(&((box_t) &(up[i])[var*box_size])[1][1][1],
                            &((box_t) &b_array[var*box_size])[xstart][ystart][z_block_size+1],
                            x_block_half, (z_block_size + 2)*(y_block_size + 2),
                            y_block_half, z_block_size + 2);
        }
    
    // DOWN -> CENTER
    if (down_f == 2)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(down[0])[var*box_size])[1][1][1],
                &((box_t) &b_array[var*box_size])[1][1][0],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                y_block_size, z_block_size + 2);
    else if (down_f == 0)
      for (int var = start; var < start + number; ++var)
        comm_eq(&((box_t) &(down[0])[var*box_size])[1][1][z_block_size],
                &((box_t) &b_array[var*box_size])[1][1][0],
                x_block_size, (z_block_size + 2)*(y_block_size + 2),
                y_block_size, z_block_size + 2);
    else if (down_f < 0) {
      down_f = -down_f - 1;
      int xstart = 1 + (down_f%2)*x_block_half;
      int ystart = 1 + (down_f/2)*y_block_half;
      for (int var = start; var < start + number; ++var)
        comm_coarsetofine(&((box_t) &(down[0])[var*box_size])[xstart][ystart][z_block_size],
                          &((box_t) &b_array[var*box_size])[1][1][0],
                          x_block_half, (z_block_size + 2)*(y_block_size + 2),
                          y_block_half, z_block_size + 2);
    }
    else
      for (int var = 0; var < start + number; ++var)
        for (int i = 0; i < 4; ++i) {
          int xstart = 1 + (i%2)*x_block_half;
          int ystart = 1 + (i/2)*y_block_half;
          comm_finetocoarse(&((box_t) &(down[i])[var*box_size])[1][1][z_block_size],
                            &((box_t) &b_array[var*box_size])[xstart][ystart][0],
                            x_block_half, (z_block_size + 2)*(y_block_size + 2),
                            y_block_half, z_block_size + 2);

          // printf("Fine to coarse %d %d\n", xstart, ystart);
        }
    

    free(north);
    free(south);
    free(east);
    free(west);
    free(up);
    free(down);
  }
}

