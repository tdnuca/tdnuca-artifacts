#include <stdio.h>
#include "stencil_kernel.h"

void printMat(int dim, int ts, double  (* restrict m)[dim][dim], int phalo, int nofirstblock) ;
void printBlock(int dim, int ts, int j, int endj, int i, int endi, double (* restrict m)[dim][dim], int phalo, int nofirstblock);

void stencil_kernel_star( int n, double (* restrict in)[n][n], double (* restrict out)[n][n], int jt, int it, int tile_size, double (* restrict weight)[2*RADIUS+1][2*RADIUS+1] ) {
   int j, i, ii, jj;
   for (j=jt; j<jt+tile_size; j++) {
      for (i=it; i<it+tile_size; i++) {
         for (jj=-RADIUS; jj<=RADIUS; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in)[j+jj][i];
         }
         for (ii=-RADIUS; ii<0; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in)[j][i+ii];
         }
         for (ii=1; ii<=RADIUS; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in)[j][i+ii];
         }
      }
   }
}

void stencil_kernel_star_3in( int n, double (* restrict in_center)[n][n],
                                     double (* restrict in_up)[n][n],
                                     double (* restrict in_down)[n][n],
                                     double (* restrict out)[n][n],
                                     int jt, int it, int tile_size, double (* restrict weight)[2*RADIUS+1][2*RADIUS+1] ) {
   int j, i, ii, jj;
   //fprintf(stderr, "jt=%d it=%d in_center %p (%p), in_up %p (%p), in_down %p (%d,%d %p)\n", jt, it, in_center, &((*in_center)[jt][it]), in_up, &((*in_up)[jt-RADIUS][it]), in_down, jt+tile_size, it, &((*in_down)[jt+tile_size][it]) );
   /* up halo */
   for (j=jt; j<jt+RADIUS; j++) {
      for (i=it; i<it+tile_size; i++) {
         /* row: center */
         for (ii=-RADIUS; ii<=RADIUS; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
         /* column, top to center, up ref */
         for (jj=-RADIUS; jj<jt-j; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_up)[j+jj][i];
         }
         /* column, top to center, center ref */
         for (jj=jt-j; jj<0; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
         /* colum, center to bottom */
         for (jj=1; jj<=RADIUS; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
      }
   }

   /* up halo */
   for (j=jt+(tile_size-RADIUS); j<jt+tile_size; j++) {
      for (i=it; i<it+tile_size; i++) {
         /* row: center */
         for (ii=-RADIUS; ii<=RADIUS; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
         /* column, center to bottom, center ref */
         for (jj=1; jj<(jt+tile_size)-j; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
         /* column, center to bottom, down ref */
         for (jj=(jt+tile_size)-j; jj<=RADIUS; jj++) {
         //fprintf(stderr, "access to down ref: base=%p, elem[%d (j=%d,jj=%d)][%d]=%p val=%lf\n", &((*in_down)[0][0]), j+jj, j, jj, i, &((*in_down)[j+jj][i]), (*in_down)[j+jj][i] );
         //fprintf(stderr, "access to down ref: base=%p, elem[%d (j=%d,jj=%d)][%d]=%p", &((*in_down)[0][0]), j+jj, j, jj, i, &((*in_down)[j+jj][i]) );
         //fprintf(stderr, " val=%lf\n", (*in_down)[j+jj][i] );
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_down)[j+jj][i];
         }
         /* colum, top to center */
         for (jj=-RADIUS; jj<0; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
      }
   }

   /* center */
   for (j=jt+RADIUS; j<(jt+tile_size)-RADIUS; j++) {
      for (i=it; i<it+tile_size; i++) {
         for (jj=-RADIUS; jj<=RADIUS; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
         for (ii=-RADIUS; ii<0; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
         for (ii=1; ii<=RADIUS; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
      }
   }
}

void stencil_kernel_star_3in_dump_in( int n, double (* restrict in_center)[n][n],
                                     double (* restrict in_up)[n][n],
                                     double (* restrict in_down)[n][n],
                                     int jt, int it, int tile_size, double (* restrict weight)[2*RADIUS+1][2*RADIUS+1] ) {
   int j, i, ii, jj;
   //fprintf(stderr, "jt=%d it=%d in_center %p (%p), in_up %p (%p), in_down %p (%d,%d %p)\n", jt, it, in_center, &((*in_center)[jt][it]), in_up, &((*in_up)[jt-RADIUS][it]), in_down, jt+tile_size, it, &((*in_down)[jt+tile_size][it]) );
   /* up halo */
   for (j=jt; j<=jt; j++) {
      for (i=it; i<=it; i++) {
         /* row: center */
         for (ii=-RADIUS; ii<=RADIUS; ii++) {
            printf("cr_ (*out)[%d][%d] += weight[%d][%d]=%lf * in_center[%d][%d]{%p}=%lf\n", j, i, ii+RADIUS, RADIUS, (*weight)[ii+RADIUS][RADIUS], j, i+ii, &(*in_center)[j][i+ii], (*in_center)[j][i+ii] );
         }
         /* column, top to center, up ref */
         for (jj=-RADIUS; jj<jt-j; jj++) {
            //(*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_up)[j+jj][i];
            printf("tc0 (*out)[%d][%d] += weight[%d][%d]=%lf * in_up[%d][%d]{%p}=%lf\n", j, i, RADIUS, jj+RADIUS, (*weight)[RADIUS][jj+RADIUS], j+jj, i, &(*in_up)[j+jj][i], (*in_up)[j+jj][i] );
         }
         /* column, top to center, center ref */
         for (jj=jt-j; jj<0; jj++) {
            //(*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
            printf("tc1 (*out)[%d][%d] += weight[%d][%d]=%lf * in_center[%d][%d]{%p}=%lf\n", j, i, RADIUS, jj+RADIUS, (*weight)[RADIUS][jj+RADIUS], j+jj, i, &(*in_center)[j+jj][i], (*in_center)[j+jj][i] );
         }
         /* colum, center to bottom */
         for (jj=1; jj<=RADIUS; jj++) {
            //(*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
            printf("bc_ (*out)[%d][%d] += weight[%d][%d]=%lf * in_center[%d][%d]{%p}=%lf\n", j, i, RADIUS, jj+RADIUS, (*weight)[RADIUS][jj+RADIUS], j+jj, i, &(*in_center)[j+jj][i], (*in_center)[j+jj][i] );
         }
      }
   }
#if 0
   /* up halo */
   for (j=jt+(tile_size-RADIUS); j<jt+tile_size; j++) {
      for (i=it; i<it+tile_size; i++) {
         /* row: center */
         for (ii=-RADIUS; ii<=RADIUS; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
         /* column, center to bottom, center ref */
         for (jj=1; jj<(jt+tile_size)-j; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
         /* column, center to bottom, down ref */
         for (jj=(jt+tile_size)-j; jj<=RADIUS; jj++) {
         //fprintf(stderr, "access to down ref: base=%p, elem[%d (j=%d,jj=%d)][%d]=%p val=%lf\n", &((*in_down)[0][0]), j+jj, j, jj, i, &((*in_down)[j+jj][i]), (*in_down)[j+jj][i] );
         //fprintf(stderr, "access to down ref: base=%p, elem[%d (j=%d,jj=%d)][%d]=%p", &((*in_down)[0][0]), j+jj, j, jj, i, &((*in_down)[j+jj][i]) );
         //fprintf(stderr, " val=%lf\n", (*in_down)[j+jj][i] );
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_down)[j+jj][i];
         }
         /* colum, top to center */
         for (jj=-RADIUS; jj<0; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
      }
   }

   /* center */
   for (j=jt+RADIUS; j<(jt+tile_size)-RADIUS; j++) {
      for (i=it; i<it+tile_size; i++) {
         for (jj=-RADIUS; jj<=RADIUS; jj++) {
            (*out)[j][i] += (*weight)[RADIUS][jj+RADIUS]*(*in_center)[j+jj][i];
         }
         for (ii=-RADIUS; ii<0; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
         for (ii=1; ii<=RADIUS; ii++) {
            (*out)[j][i] += (*weight)[ii+RADIUS][RADIUS]*(*in_center)[j][i+ii];
         }
      }
   }
   #endif
   printf("oh %d %d, %d %d\n", jt, jt+tile_size, it, it+tile_size);
   printBlock(n, tile_size, jt, jt+tile_size, it, it+tile_size, in_center, 0, 0);
}
