// compute_max changed; now not indexed through ibool
// NDIM moved to rightmost dimension
/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

//
// All the arrays below use static memory allocation,
// using constant sizes defined in values_from_mesher.h.
// This is done purposely to improve performance (Fortran compilers
// can optimize much more when the size of the loops and arrays
// is known at compile time).
// NGLLX, NGLLY and NGLLZ are set equal to 5,
// therefore each element contains NGLLX * NGLLY * NGLLZ = 125 points.
//

//
// All the calculations are done in single precision.
// We do not need double precision in SPECFEM3D.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>



// use the Deville et al. (2002) inlined products or not
#define USE_DEVILLE_INLINED_PRODUCTS

// include values created by the mesher
// done for performance only using static allocation to allow for loop unrolling
#include "./values_from_mesher_C.h"

// constant value of the time step in the main time loop
#define deltatover2 0.5f*deltat
#define deltatsqover2 0.5f*deltat*deltat

// for the source time function
#define pi 3.141592653589793f
#define f0 (1.f / 50.f)
#define t0 (1.2f / f0)
#define a pi*pi*f0*f0

// number of GLL integration points in each direction of an element (degree plus one)
//using enums to get real C constants
//because defines cannot be used for OmpSs pragmas
enum { NGLLX = 5 };
enum { NGLLY = 5 };
enum { NGLLZ = 5 };

// for the Deville et al. (2002) inlined products
#define NGLL2  25 // NGLLX^2

// 3-D simulation
//using enums to get real C constants
//because defines cannot be used for OmpSs pragmas
enum { NDIM = 3 };

// displacement threshold above which we consider that the code became unstable
#define STABILITY_THRESHOLD 1.e+25f

//#define NTSTEP_BETWEEN_OUTPUT_INFO 1000
#define NTSTEP_BETWEEN_OUTPUT_INFO 100

// approximate density of the geophysical medium in which the source is located
// this value is only a constant scaling factor therefore it does not really matter
#define rho 4500.f

// CellSs or SMPSs: Increase block size in order to reduce overhead.
// DK DK we should see if there are roundoff problems here when NSPEC is not a multiple of NUMBER_OF_THREADS_USED
//#define NUMBER_OF_THREADS_USED 8
//#define BS_NGLOB (NGLOB / NUMBER_OF_THREADS_USED)
//#define BS_NSPEC (NSPEC / NUMBER_OF_THREADS_USED)
// other values from Rosa
//#define BS_NGLOB 10000
//#define BS_NSPEC 200
//#define PROP_THREADS  //No se perque amb aixo no em funciona :-( 
//#define BS_NGLOB 4000
//#define BS_NGLOB 4614
#define BS_NGLOB 14765 //vdimic
#define BS_NSPEC 3400 //vdimic
//#define BS_NSPEC 30
//#define BS_NGLOB 2325
//#define BS_NSPEC 10
//#define BS_NGLOB 17
//Tricky, to obtaine 4 chuncks per thread with 16 procs
//#define BS_NSPEC 34

// define some constants to group the arrays to reduce the number of arguments
// to send to the tasks
#define XIX 0
#define XIY 1
#define XIZ 2
#define ETAX 3
#define ETAY 4
#define ETAZ 5
#define GAMMAX 6
#define GAMMAY 7
#define GAMMAZ 8

#define KAPPA 0
#define MU 1

#define X 0
#define Y 1
#define Z 2

#define YZ 0
#define XZ 1
#define XY 2

#define FLAG_hprime_xx 0
#define FLAG_hprime_xxT 1
#define FLAG_hprimewgll_xx 2
#define FLAG_hprimewgll_xxT 3
#define FLAG_wgllwgll_YZ 4
#define FLAG_wgllwgll_XZ 5
#define FLAG_wgllwgll_XY 6



long t_start,t_end;

// declare all the functions
static long usecs ();

static void clear (int actual_size, float displ[actual_size][NDIM], float veloc[actual_size][NDIM], float accel[actual_size][NDIM]);

/*static void gather (int actual_size,
             float displ[NGLOB][NDIM],
             int ibool[actual_size][NGLLZ][NGLLY][NGLLX],
             float dummy_loc[actual_size][NGLLZ][NGLLY][NGLLX][NDIM]);

static void scatter (int actual_size,
              int ibool[actual_size][NGLLZ][NGLLY][NGLLX],
              float sum_terms[actual_size][NGLLZ][NGLLY][NGLLX][NDIM],
              float accel[NGLOB][NDIM]);

static void process_element(
  int actual_size,
  float dummy_loc[actual_size][NGLLZ][NGLLY][NGLLX][NDIM],
#ifdef USE_DEVILLE_INLINED_PRODUCTS
  float all_matrices[NDIM+4][NGLLX][NGLLX],
#else
  float hprime_xx[NGLLX][NGLLX],
  float hprimewgll_xx[NGLLX][NGLLX],
  float wgllwgll_all[NGLLZ][NGLLY][NDIM],
#endif
  float jacobian_matrix[actual_size][NGLLZ][NGLLY][NGLLX][NDIM*NDIM],
  float kappa_and_mu[actual_size][NGLLZ][NGLLY][NGLLX][2],
  float sum_terms[actual_size][NGLLZ][NGLLY][NGLLX][NDIM]);

static void update_disp_vel (int actual_size,
                      float displ[actual_size][NDIM],
                      float veloc[actual_size][NDIM],
                      float accel[actual_size][NDIM]);

static void update_acc_vel (int actual_size,
                     float accel[actual_size][NDIM],
                     float veloc[actual_size][NDIM],
                     float rmass_inverse[actual_size]);

static void compute_max(int actual_size,
                 float displ[NGLOB][NDIM],
//                 int ibool[actual_size][NGLLZ][NGLLY][NGLLX],
                 float *max);

*/
////////////////////////////////////////////////////////////////////////
//                                   TASKS                            //
////////////////////////////////////////////////////////////////////////


#pragma omp task out(displ[0;actual_size][0;NDIM], veloc[0;actual_size][0;NDIM], accel[0;actual_size][0;NDIM])
void clear (int actual_size, float displ[actual_size][NDIM], float veloc[actual_size][NDIM], float accel[actual_size][NDIM])
{
  int i;
  for (i=0;i<actual_size;i++) {
    displ[i][X] = 0.f;
    displ[i][Y] = 0.f;
    displ[i][Z] = 0.f;

    veloc[i][X] = 0.f;
    veloc[i][Y] = 0.f;
    veloc[i][Z] = 0.f;

    accel[i][X] = 0.f;
    accel[i][Y] = 0.f;
    accel[i][Z] = 0.f;
  }

}


////////////////////////////////////////////////////////////////
// Updates
///////////////////////////////////////////////////////////////

// update the displacement and velocity vectors and clear the acceleration
// vector for the assinged chunck

#pragma omp task inout(displ[0;actual_size][0;NDIM], veloc[0;actual_size][0;NDIM], accel[0;actual_size][0;NDIM])
void update_disp_vel (int actual_size,
                      float displ[actual_size][NDIM],
                      float veloc[actual_size][NDIM],
                      float accel[actual_size][NDIM])
{
   int i;

//printf("update_disp_vel displ=%d, veloc=%d, accel=%d\n", NDIM*actual_size*sizeof(float), NDIM*actual_size*sizeof(float), NDIM*actual_size*sizeof(float));
// DK DK we CANNOT define this with the [NDIM] index first, as would be natural in C
// DK DK in order to have the fastest index [i] on the right, because these arrays
// DK DK are defined as accel[actual_size][NDIM] in the tasks but as accel[NGLOB][NDIM]
// DK DK in other routines and therefore the memory chunks would not correspond.
// DK DK In other words, in any array with a dimension [actual_size] that dimension must
// DK DK always be the leftmost index.
   for (i=0; i<actual_size; i++){
    displ[i][X] += deltat*veloc[i][X] + deltatsqover2*accel[i][X];
    displ[i][Y] += deltat*veloc[i][Y] + deltatsqover2*accel[i][Y];
    displ[i][Z] += deltat*veloc[i][Z] + deltatsqover2*accel[i][Z];

    veloc[i][X] += deltatover2*accel[i][X];
    veloc[i][Y] += deltatover2*accel[i][Y];
    veloc[i][Z] += deltatover2*accel[i][Z];

    accel[i][X] = 0.f;
    accel[i][Y] = 0.f;
    accel[i][Z] = 0.f;
  }
}

// update the acceleration and velocity vectors on assigned chunk of points

#pragma omp task in(rmass_inverse[0;actual_size]) inout(accel[0;actual_size][0;NDIM], veloc[0;actual_size][0;NDIM])
void update_acc_vel (int actual_size,
                     float accel[actual_size][NDIM],
                     float veloc[actual_size][NDIM],
                     float rmass_inverse[actual_size],
		     int iglob_s,
		     float time)
{
   int i;

  if (iglob_s != -1) accel[iglob_s][Z] += 1.e4f * (1.f - 2.f*a*(time-t0)*(time-t0)) * expf(-a*(time-t0)*(time-t0)) / (rho * rmass_inverse[iglob_s]);

   for (i=0; i<actual_size; i++){
    accel[i][X] *= rmass_inverse[i];
    accel[i][Y] *= rmass_inverse[i];
    accel[i][Z] *= rmass_inverse[i];

    veloc[i][X] += deltatover2*accel[i][X];
    veloc[i][Y] += deltatover2*accel[i][Y];
    veloc[i][Z] += deltatover2*accel[i][Z];
  }
}

/// playing around ///


void zones (int connectivity[NSPEC/BS_NSPEC+1][NGLOB/BS_NGLOB+1],int ibool[NSPEC][NGLLZ][NGLLY][NGLLX]){
int zona, new_zona;
//int accesses;
int elem, i, j, k, l;
int cacho;

//int accesos[NGLOB/BS_NGLOB+1];
//zona = ibool[0][0][0][0]/BS_NGLOB;
//printf("Zone %d:", zona);
//accesses=1;
for (j=0; j<NSPEC/BS_NSPEC+1;j++)
   for (i=0; i<NGLOB/BS_NGLOB+1; i++) 
	connectivity[j][i]=0;
cacho=0;
for (elem=0; elem<NSPEC; elem++){ 
	if ((elem % BS_NSPEC ==0)&& (elem>0)){
	 cacho ++;
//	 printf ("\n");
//	 printf ("------------------------------\n");
//	 for (l=0; l < NGLOB/BS_NGLOB; l++) {
//		if (accesos[l]!=0) printf("zona %d accesos %d \n", l, accesos[l]);
//		if (accesos[l]!=0) printf("X");
//		else printf(" ");
//		accesos[l]=0;
//	 }
	}
      for (k=0; k<NGLLZ;k++)
        for (j=0;j<NGLLY;j++)
          for (i=0;i<NGLLX;i++) {
              new_zona = ibool[elem][k][j][i]/BS_NGLOB;
	      connectivity[cacho][new_zona]++;
//		printf("ibool %d \n", ibool[elem][k][j][i]);
/*	      if (zona == new_zona) accesos[new_zona]++;
	      else {
		zona = new_zona;
		printf(" %d accesses\n", accesses);
		accesses=1;
		printf("Zone %d:", zona);
	      }*/
	//	 printf("new_zona = %d \n"); 
	}	
}

for (i=0; i < NSPEC/BS_NSPEC+1; i++){
	 for (l=0; l < NGLOB/BS_NGLOB; l++) {
		if (connectivity[i][l]!=0) printf("X");
		else printf(" ");
	}
	printf("\n");
}
	  
}

/////////////////////////////////////////////////////////////////////
//   gather - scatter
/////////////////////////////////////////////////////////////////////


// localize data for the element from the global vectors to the local mesh
// using indirect addressing (contained in array ibool)

//#pragma omp task input(elem_size, displ, ibool, node_size, zona) inout(dummy_loc) 
#pragma omp task in(displ[0;node_size], ibool[0;elem_size][0;NGLLZ][0;NGLLY][0;NGLLX]) out(dummy_loc[0;elem_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM])
void gather (int elem_size,
	     int node_size,
             float displ[node_size][NDIM],
             int ibool[elem_size][NGLLZ][NGLLY][NGLLX],
             float dummy_loc[elem_size][NGLLZ][NGLLY][NGLLX][NDIM],
	     int zona)
{
   int i,j,k,iglob,elem;
   int zona_local;

   for (elem=0; elem<elem_size; elem++) {
      for (k=0; k<NGLLZ;k++)
        for (j=0;j<NGLLY;j++)
          for (i=0;i<NGLLX;i++) {
	      zona_local = ibool[elem][k][j][i]/BS_NGLOB;
	      if (zona_local == zona) {
//		printf("elem k j i %d %d %d %d \n", elem, k, j, i);
                iglob = ibool[elem][k][j][i]%BS_NGLOB;
              	dummy_loc[elem][k][j][i][X] = displ[iglob][X];
              	dummy_loc[elem][k][j][i][Y] = displ[iglob][Y];
              	dummy_loc[elem][k][j][i][Z] = displ[iglob][Z];
	      }
          }
   }
}

// scattered update must be atomic because we did not use mesh coloring to make
// mesh subsets independent. We could use mesh coloring instead, as in Figure 6 of
// Dimitri Komatitsch, David Michea and Gordon Erlebacher,
// Porting a high-order finite-element earthquake modeling application to NVIDIA graphics cards using CUDA,
// Journal of Parallel and Distributed Computing,
// vol. 69(5), p. 451-460, doi: 10.1016/j.jpdc.2009.01.006 (2009).
// http://www.univ-pau.fr/~dkomati1/published_papers/GPGPU_JPDC_2009.pdf

// sum contributions from the element to the global mesh using indirect addressing
// Sequentiality imposed by dependence on whole accel

//#pragma omp task input(elem_size, sum_terms,ibool, node_size, zona) inout(accel)
#pragma omp task in(sum_terms[0;elem_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM], ibool[0;elem_size][0;NGLLZ][0;NGLLY][0;NGLLX]) out(accel[0;node_size][0;NDIM])
void scatter (int elem_size,
	      int node_size,
              int ibool[elem_size][NGLLZ][NGLLY][NGLLX],
              float sum_terms[elem_size][NGLLZ][NGLLY][NGLLX][NDIM],
              float accel[node_size][NDIM],
	      int zona)
{
   int i,j,k,iglob, elem;
   int zona_local;
   
   for (elem=0; elem<elem_size; elem++) {
     for (k=0;k<NGLLZ;k++) {
        for (j=0;j<NGLLY;j++) {
           for (i=0;i<NGLLX;i++) {
             zona_local = ibool[elem][k][j][i]/BS_NGLOB;
	     if (zona_local == zona) {
		iglob = ibool[elem][k][j][i]%BS_NGLOB;
#pragma mcc verbatim start 
#pragma omp atomic
#pragma mcc verbatim end
             accel[iglob][X] += sum_terms[elem][k][j][i][X];
#pragma mcc verbatim start 
#pragma omp atomic
#pragma mcc verbatim end
             accel[iglob][Y] += sum_terms[elem][k][j][i][Y];
#pragma mcc verbatim start 
#pragma omp atomic
#pragma mcc verbatim end
             accel[iglob][Z] += sum_terms[elem][k][j][i][Z];
  	      }
           }
         }
      }

   }
}

////////////////////////////////////////////////////////////////////
//      Element
///////////////////////////////////////////////////////////////////

#ifdef USE_DEVILLE_INLINED_PRODUCTS
#pragma omp task in (dummy_loc[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM], all_matrices[0;NDIM+4][0;NGLLX][0;NGLLX], jacobian_matrix[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;2]) \
                in (kappa_and_mu[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;2]) out (sum_terms[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM])
void process_element(
  int actual_size,
  float dummy_loc[actual_size][NGLLZ][NGLLY][NGLLX][NDIM],
  float all_matrices[NDIM+4][NGLLX][NGLLX],
  float jacobian_matrix[actual_size][NGLLZ][NGLLY][NGLLX][NDIM*NDIM],
  float kappa_and_mu[actual_size][NGLLZ][NGLLY][NGLLX][2],
  float sum_terms[actual_size][NGLLZ][NGLLY][NGLLX][NDIM]) {

  float tempx1[NGLLZ][NGLLY][NGLLX];
  float tempx2[NGLLZ][NGLLY][NGLLX];
  float tempx3[NGLLZ][NGLLY][NGLLX];
  float tempy1[NGLLZ][NGLLY][NGLLX];
  float tempy2[NGLLZ][NGLLY][NGLLX];
  float tempy3[NGLLZ][NGLLY][NGLLX];
  float tempz1[NGLLZ][NGLLY][NGLLX];
  float tempz2[NGLLZ][NGLLY][NGLLX];
  float tempz3[NGLLZ][NGLLY][NGLLX];

  float newtempx1[NGLLZ][NGLLY][NGLLX];
  float newtempx2[NGLLZ][NGLLY][NGLLX];
  float newtempx3[NGLLZ][NGLLY][NGLLX];
  float newtempy1[NGLLZ][NGLLY][NGLLX];
  float newtempy2[NGLLZ][NGLLY][NGLLX];
  float newtempy3[NGLLZ][NGLLY][NGLLX];
  float newtempz1[NGLLZ][NGLLY][NGLLX];
  float newtempz2[NGLLZ][NGLLY][NGLLX];
  float newtempz3[NGLLZ][NGLLY][NGLLX];

  int i, j, k, elem;
  float xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  float duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  float duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  float duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  float sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  float fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;

// DK DK pointers to simulate a "union" statement by pointing to the same memory block
  float *tempx1_2D_25_5, *tempy1_2D_25_5, *tempz1_2D_25_5;
  float *tempx3_2D_5_25, *tempy3_2D_5_25, *tempz3_2D_5_25;
  float *newtempx1_2D_25_5, *newtempy1_2D_25_5, *newtempz1_2D_25_5;
  float *newtempx3_2D_5_25, *newtempy3_2D_5_25, *newtempz3_2D_5_25;
  float *dummyx_loc_2D, *dummyy_loc_2D, *dummyz_loc_2D;

//printf("procces_element actual_size = %d\n", actual_size);

// DK DK mapping for 2D arrays of size [25][5] or [5][25] to an offset in a 1D linear memory block
#define map_25_5(i,j)   ((j) +  5*(i))
#define map_5_25(i,j)   ((j) + 25*(i))

// DK DK assign the pointers to the common memory block
  tempx1_2D_25_5 = &tempx1[0][0][0];
  tempy1_2D_25_5 = &tempy1[0][0][0];
  tempz1_2D_25_5 = &tempz1[0][0][0];

  tempx3_2D_5_25 = &tempx3[0][0][0];
  tempy3_2D_5_25 = &tempy3[0][0][0];
  tempz3_2D_5_25 = &tempz3[0][0][0];

  newtempx1_2D_25_5 = &newtempx1[0][0][0];
  newtempy1_2D_25_5 = &newtempy1[0][0][0];
  newtempz1_2D_25_5 = &newtempz1[0][0][0];

  newtempx3_2D_5_25 = &newtempx3[0][0][0];
  newtempy3_2D_5_25 = &newtempy3[0][0][0];
  newtempz3_2D_5_25 = &newtempz3[0][0][0];

  for (elem=0; elem<actual_size; elem++) {

// DK DK assign the pointers to the common memory block
  dummyx_loc_2D = &dummy_loc[elem][0][0][0][X];
  dummyy_loc_2D = &dummy_loc[elem][0][0][0][Y];
  dummyz_loc_2D = &dummy_loc[elem][0][0][0][Z];

// subroutines adapted from Deville, Fischer and Mund, High-order methods
// for incompressible fluid flow, Cambridge University Press (2002),
// pages 386 and 389 and Figure 8.3.1
  for (i=0;i<NGLLX;i++) {
    for (j=0;j<NGLL2;j++) {
      tempx1_2D_25_5[map_25_5(j,i)] = all_matrices[FLAG_hprime_xx][0][i]*dummyx_loc_2D[map_25_5(j,0)*3] +
                                      all_matrices[FLAG_hprime_xx][1][i]*dummyx_loc_2D[map_25_5(j,1)*3] +
                                      all_matrices[FLAG_hprime_xx][2][i]*dummyx_loc_2D[map_25_5(j,2)*3] +
                                      all_matrices[FLAG_hprime_xx][3][i]*dummyx_loc_2D[map_25_5(j,3)*3] +
                                      all_matrices[FLAG_hprime_xx][4][i]*dummyx_loc_2D[map_25_5(j,4)*3];

      tempy1_2D_25_5[map_25_5(j,i)] = all_matrices[FLAG_hprime_xx][0][i]*dummyy_loc_2D[map_25_5(j,0)*3] +
                                      all_matrices[FLAG_hprime_xx][1][i]*dummyy_loc_2D[map_25_5(j,1)*3] +
                                      all_matrices[FLAG_hprime_xx][2][i]*dummyy_loc_2D[map_25_5(j,2)*3] +
                                      all_matrices[FLAG_hprime_xx][3][i]*dummyy_loc_2D[map_25_5(j,3)*3] +
                                      all_matrices[FLAG_hprime_xx][4][i]*dummyy_loc_2D[map_25_5(j,4)*3];

      tempz1_2D_25_5[map_25_5(j,i)] = all_matrices[FLAG_hprime_xx][0][i]*dummyz_loc_2D[map_25_5(j,0)*3] +
                                      all_matrices[FLAG_hprime_xx][1][i]*dummyz_loc_2D[map_25_5(j,1)*3] +
                                      all_matrices[FLAG_hprime_xx][2][i]*dummyz_loc_2D[map_25_5(j,2)*3] +
                                      all_matrices[FLAG_hprime_xx][3][i]*dummyz_loc_2D[map_25_5(j,3)*3] +
                                      all_matrices[FLAG_hprime_xx][4][i]*dummyz_loc_2D[map_25_5(j,4)*3];
    }
  }

  for (k=0;k<NGLLZ;k++) {
    for (j=0;j<NGLLX;j++) {
      for (i=0;i<NGLLX;i++) {

        tempx2[k][j][i] = dummy_loc[elem][k][0][i][X]*all_matrices[FLAG_hprime_xxT][j][0] +
                          dummy_loc[elem][k][1][i][X]*all_matrices[FLAG_hprime_xxT][j][1] +
                          dummy_loc[elem][k][2][i][X]*all_matrices[FLAG_hprime_xxT][j][2] +
                          dummy_loc[elem][k][3][i][X]*all_matrices[FLAG_hprime_xxT][j][3] +
                          dummy_loc[elem][k][4][i][X]*all_matrices[FLAG_hprime_xxT][j][4];

        tempy2[k][j][i] = dummy_loc[elem][k][0][i][Y]*all_matrices[FLAG_hprime_xxT][j][0] +
                          dummy_loc[elem][k][1][i][Y]*all_matrices[FLAG_hprime_xxT][j][1] +
                          dummy_loc[elem][k][2][i][Y]*all_matrices[FLAG_hprime_xxT][j][2] +
                          dummy_loc[elem][k][3][i][Y]*all_matrices[FLAG_hprime_xxT][j][3] +
                          dummy_loc[elem][k][4][i][Y]*all_matrices[FLAG_hprime_xxT][j][4];

        tempz2[k][j][i] = dummy_loc[elem][k][0][i][Z]*all_matrices[FLAG_hprime_xxT][j][0] +
                          dummy_loc[elem][k][1][i][Z]*all_matrices[FLAG_hprime_xxT][j][1] +
                          dummy_loc[elem][k][2][i][Z]*all_matrices[FLAG_hprime_xxT][j][2] +
                          dummy_loc[elem][k][3][i][Z]*all_matrices[FLAG_hprime_xxT][j][3] +
                          dummy_loc[elem][k][4][i][Z]*all_matrices[FLAG_hprime_xxT][j][4];

        }
      }
    }

  for (j=0;j<NGLLX;j++) {
    for (i=0;i<NGLL2;i++) {
      tempx3_2D_5_25[map_5_25(j,i)] = dummyx_loc_2D[map_5_25(0,i)*3]*all_matrices[FLAG_hprime_xxT][j][0] +
                                      dummyx_loc_2D[map_5_25(1,i)*3]*all_matrices[FLAG_hprime_xxT][j][1] +
                                      dummyx_loc_2D[map_5_25(2,i)*3]*all_matrices[FLAG_hprime_xxT][j][2] +
                                      dummyx_loc_2D[map_5_25(3,i)*3]*all_matrices[FLAG_hprime_xxT][j][3] +
                                      dummyx_loc_2D[map_5_25(4,i)*3]*all_matrices[FLAG_hprime_xxT][j][4];

      tempy3_2D_5_25[map_5_25(j,i)] = dummyy_loc_2D[map_5_25(0,i)*3]*all_matrices[FLAG_hprime_xxT][j][0] +
                                      dummyy_loc_2D[map_5_25(1,i)*3]*all_matrices[FLAG_hprime_xxT][j][1] +
                                      dummyy_loc_2D[map_5_25(2,i)*3]*all_matrices[FLAG_hprime_xxT][j][2] +
                                      dummyy_loc_2D[map_5_25(3,i)*3]*all_matrices[FLAG_hprime_xxT][j][3] +
                                      dummyy_loc_2D[map_5_25(4,i)*3]*all_matrices[FLAG_hprime_xxT][j][4];

      tempz3_2D_5_25[map_5_25(j,i)] = dummyz_loc_2D[map_5_25(0,i)*3]*all_matrices[FLAG_hprime_xxT][j][0] +
                                      dummyz_loc_2D[map_5_25(1,i)*3]*all_matrices[FLAG_hprime_xxT][j][1] +
                                      dummyz_loc_2D[map_5_25(2,i)*3]*all_matrices[FLAG_hprime_xxT][j][2] +
                                      dummyz_loc_2D[map_5_25(3,i)*3]*all_matrices[FLAG_hprime_xxT][j][3] +
                                      dummyz_loc_2D[map_5_25(4,i)*3]*all_matrices[FLAG_hprime_xxT][j][4];

    }
  }

  for (k=0;k<NGLLZ;k++) {
    for (j=0;j<NGLLX;j++) {
      for (i=0;i<NGLLX;i++) {

// compute derivatives of ux, uy and uz with respect to x, y and z
          xixl = jacobian_matrix[elem][k][j][i][XIX];
          xiyl = jacobian_matrix[elem][k][j][i][XIY];
          xizl = jacobian_matrix[elem][k][j][i][XIZ];
          etaxl = jacobian_matrix[elem][k][j][i][ETAX];
          etayl = jacobian_matrix[elem][k][j][i][ETAY];
          etazl = jacobian_matrix[elem][k][j][i][ETAZ];
          gammaxl = jacobian_matrix[elem][k][j][i][GAMMAX];
          gammayl = jacobian_matrix[elem][k][j][i][GAMMAY];
          gammazl = jacobian_matrix[elem][k][j][i][GAMMAZ];
          jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)-xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl));

          duxdxl = xixl*tempx1[k][j][i] + etaxl*tempx2[k][j][i] + gammaxl*tempx3[k][j][i];
          duxdyl = xiyl*tempx1[k][j][i] + etayl*tempx2[k][j][i] + gammayl*tempx3[k][j][i];
          duxdzl = xizl*tempx1[k][j][i] + etazl*tempx2[k][j][i] + gammazl*tempx3[k][j][i];

          duydxl = xixl*tempy1[k][j][i] + etaxl*tempy2[k][j][i] + gammaxl*tempy3[k][j][i];
          duydyl = xiyl*tempy1[k][j][i] + etayl*tempy2[k][j][i] + gammayl*tempy3[k][j][i];
          duydzl = xizl*tempy1[k][j][i] + etazl*tempy2[k][j][i] + gammazl*tempy3[k][j][i];

          duzdxl = xixl*tempz1[k][j][i] + etaxl*tempz2[k][j][i] + gammaxl*tempz3[k][j][i];
          duzdyl = xiyl*tempz1[k][j][i] + etayl*tempz2[k][j][i] + gammayl*tempz3[k][j][i];
          duzdzl = xizl*tempz1[k][j][i] + etazl*tempz2[k][j][i] + gammazl*tempz3[k][j][i];

// precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl;
          duxdxl_plus_duzdzl = duxdxl + duzdzl;
          duydyl_plus_duzdzl = duydyl + duzdzl;
          duxdyl_plus_duydxl = duxdyl + duydxl;
          duzdxl_plus_duxdzl = duzdxl + duxdzl;
          duzdyl_plus_duydzl = duzdyl + duydzl;

// compute isotropic elements
          kappal = kappa_and_mu[elem][k][j][i][KAPPA];
          mul = kappa_and_mu[elem][k][j][i][MU];

          lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
          lambdal = lambdalplus2mul - 2.f*mul;

// compute stress sigma
          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

          sigma_xy = mul*duxdyl_plus_duydxl;
          sigma_xz = mul*duzdxl_plus_duxdzl;
          sigma_yz = mul*duzdyl_plus_duydzl;

// form dot product with test vector
          tempx1[k][j][i] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl);
          tempy1[k][j][i] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl);
          tempz1[k][j][i] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

          tempx2[k][j][i] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl);
          tempy2[k][j][i] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl);
          tempz2[k][j][i] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

          tempx3[k][j][i] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl);
          tempy3[k][j][i] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl);
          tempz3[k][j][i] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);

        }
      }
    }

  for (j=0;j<NGLL2;j++) {
    for (i=0;i<NGLLX;i++) {
      newtempx1_2D_25_5[map_25_5(j,i)] = all_matrices[FLAG_hprimewgll_xxT][0][i]*tempx1_2D_25_5[map_25_5(j,0)] +
                                         all_matrices[FLAG_hprimewgll_xxT][1][i]*tempx1_2D_25_5[map_25_5(j,1)] +
                                         all_matrices[FLAG_hprimewgll_xxT][2][i]*tempx1_2D_25_5[map_25_5(j,2)] +
                                         all_matrices[FLAG_hprimewgll_xxT][3][i]*tempx1_2D_25_5[map_25_5(j,3)] +
                                         all_matrices[FLAG_hprimewgll_xxT][4][i]*tempx1_2D_25_5[map_25_5(j,4)];

      newtempy1_2D_25_5[map_25_5(j,i)] = all_matrices[FLAG_hprimewgll_xxT][0][i]*tempy1_2D_25_5[map_25_5(j,0)] +
                                         all_matrices[FLAG_hprimewgll_xxT][1][i]*tempy1_2D_25_5[map_25_5(j,1)] +
                                         all_matrices[FLAG_hprimewgll_xxT][2][i]*tempy1_2D_25_5[map_25_5(j,2)] +
                                         all_matrices[FLAG_hprimewgll_xxT][3][i]*tempy1_2D_25_5[map_25_5(j,3)] +
                                         all_matrices[FLAG_hprimewgll_xxT][4][i]*tempy1_2D_25_5[map_25_5(j,4)];

      newtempz1_2D_25_5[map_25_5(j,i)] = all_matrices[FLAG_hprimewgll_xxT][0][i]*tempz1_2D_25_5[map_25_5(j,0)] +
                                         all_matrices[FLAG_hprimewgll_xxT][1][i]*tempz1_2D_25_5[map_25_5(j,1)] +
                                         all_matrices[FLAG_hprimewgll_xxT][2][i]*tempz1_2D_25_5[map_25_5(j,2)] +
                                         all_matrices[FLAG_hprimewgll_xxT][3][i]*tempz1_2D_25_5[map_25_5(j,3)] +
                                         all_matrices[FLAG_hprimewgll_xxT][4][i]*tempz1_2D_25_5[map_25_5(j,4)];
    }
  }

  for (k=0;k<NGLLZ;k++) {
    for (j=0;j<NGLLX;j++) {
      for (i=0;i<NGLLX;i++) {
        newtempx2[k][j][i] = tempx2[k][0][i]*all_matrices[FLAG_hprimewgll_xx][j][0] +
                             tempx2[k][1][i]*all_matrices[FLAG_hprimewgll_xx][j][1] +
                             tempx2[k][2][i]*all_matrices[FLAG_hprimewgll_xx][j][2] +
                             tempx2[k][3][i]*all_matrices[FLAG_hprimewgll_xx][j][3] +
                             tempx2[k][4][i]*all_matrices[FLAG_hprimewgll_xx][j][4];

        newtempy2[k][j][i] = tempy2[k][0][i]*all_matrices[FLAG_hprimewgll_xx][j][0] +
                             tempy2[k][1][i]*all_matrices[FLAG_hprimewgll_xx][j][1] +
                             tempy2[k][2][i]*all_matrices[FLAG_hprimewgll_xx][j][2] +
                             tempy2[k][3][i]*all_matrices[FLAG_hprimewgll_xx][j][3] +
                             tempy2[k][4][i]*all_matrices[FLAG_hprimewgll_xx][j][4];

        newtempz2[k][j][i] = tempz2[k][0][i]*all_matrices[FLAG_hprimewgll_xx][j][0] +
                             tempz2[k][1][i]*all_matrices[FLAG_hprimewgll_xx][j][1] +
                             tempz2[k][2][i]*all_matrices[FLAG_hprimewgll_xx][j][2] +
                             tempz2[k][3][i]*all_matrices[FLAG_hprimewgll_xx][j][3] +
                             tempz2[k][4][i]*all_matrices[FLAG_hprimewgll_xx][j][4];
      }
    }
  }

  for (j=0;j<NGLLX;j++) {
    for (i=0;i<NGLL2;i++) {
      newtempx3_2D_5_25[map_5_25(j,i)] = tempx3_2D_5_25[map_5_25(0,i)]*all_matrices[FLAG_hprimewgll_xx][j][0] +
                                         tempx3_2D_5_25[map_5_25(1,i)]*all_matrices[FLAG_hprimewgll_xx][j][1] +
                                         tempx3_2D_5_25[map_5_25(2,i)]*all_matrices[FLAG_hprimewgll_xx][j][2] +
                                         tempx3_2D_5_25[map_5_25(3,i)]*all_matrices[FLAG_hprimewgll_xx][j][3] +
                                         tempx3_2D_5_25[map_5_25(4,i)]*all_matrices[FLAG_hprimewgll_xx][j][4];

      newtempy3_2D_5_25[map_5_25(j,i)] = tempy3_2D_5_25[map_5_25(0,i)]*all_matrices[FLAG_hprimewgll_xx][j][0] +
                                         tempy3_2D_5_25[map_5_25(1,i)]*all_matrices[FLAG_hprimewgll_xx][j][1] +
                                         tempy3_2D_5_25[map_5_25(2,i)]*all_matrices[FLAG_hprimewgll_xx][j][2] +
                                         tempy3_2D_5_25[map_5_25(3,i)]*all_matrices[FLAG_hprimewgll_xx][j][3] +
                                         tempy3_2D_5_25[map_5_25(4,i)]*all_matrices[FLAG_hprimewgll_xx][j][4];

      newtempz3_2D_5_25[map_5_25(j,i)] = tempz3_2D_5_25[map_5_25(0,i)]*all_matrices[FLAG_hprimewgll_xx][j][0] +
                                         tempz3_2D_5_25[map_5_25(1,i)]*all_matrices[FLAG_hprimewgll_xx][j][1] +
                                         tempz3_2D_5_25[map_5_25(2,i)]*all_matrices[FLAG_hprimewgll_xx][j][2] +
                                         tempz3_2D_5_25[map_5_25(3,i)]*all_matrices[FLAG_hprimewgll_xx][j][3] +
                                         tempz3_2D_5_25[map_5_25(4,i)]*all_matrices[FLAG_hprimewgll_xx][j][4];
    }
  }

    for (k=0;k<NGLLZ;k++) {
      for (j=0;j<NGLLY;j++) {
        for (i=0;i<NGLLX;i++) {

          fac1 = all_matrices[FLAG_wgllwgll_YZ][k][j];
          fac2 = all_matrices[FLAG_wgllwgll_XZ][k][i];
          fac3 = all_matrices[FLAG_wgllwgll_XY][j][i];

          sum_terms[elem][k][j][i][X] = - (fac1*newtempx1[k][j][i] + fac2*newtempx2[k][j][i] + fac3*newtempx3[k][j][i]);
          sum_terms[elem][k][j][i][Y] = - (fac1*newtempy1[k][j][i] + fac2*newtempy2[k][j][i] + fac3*newtempy3[k][j][i]);
          sum_terms[elem][k][j][i][Z] = - (fac1*newtempz1[k][j][i] + fac2*newtempz2[k][j][i] + fac3*newtempz3[k][j][i]);

        }
      }
    }

  }

}
#else // of USE_DEVILLE_INLINED_PRODUCTS
#pragma omp task in (dummy_loc[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM], hprimewgll_xx[0;NGLLX][0;NGLLX]hprime_xx, [0;NGLLX][0;NGLLX], wgllwgll_all[0;NGLLZ][0;NGLLY][0;NDIM]) \
                 in (jacobian_matrix[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM*NDIM], kappa_and_mu[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;2]) \
                 out (sum_terms[0;actual_size][0;NGLLZ][0;NGLLY][0;NGLLX][0;NDIM])
void process_element(
  int actual_size,
  float dummy_loc[actual_size][NGLLZ][NGLLY][NGLLX][NDIM],
  float hprime_xx[NGLLX][NGLLX],
  float hprimewgll_xx[NGLLX][NGLLX],
  float wgllwgll_all[NGLLZ][NGLLY][NDIM],
  float jacobian_matrix[actual_size][NGLLZ][NGLLY][NGLLX][NDIM*NDIM],
  float kappa_and_mu[actual_size][NGLLZ][NGLLY][NGLLX][2],
  float sum_terms[actual_size][NGLLZ][NGLLY][NGLLX][NDIM]) {

  float tempx1[NGLLZ][NGLLY][NGLLX];
  float tempx2[NGLLZ][NGLLY][NGLLX];
  float tempx3[NGLLZ][NGLLY][NGLLX];
  float tempy1[NGLLZ][NGLLY][NGLLX];
  float tempy2[NGLLZ][NGLLY][NGLLX];
  float tempy3[NGLLZ][NGLLY][NGLLX];
  float tempz1[NGLLZ][NGLLY][NGLLX];
  float tempz2[NGLLZ][NGLLY][NGLLX];
  float tempz3[NGLLZ][NGLLY][NGLLX];

  int i, j, k, l, elem;
  float xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  float duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  float duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  float duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  float sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  float hp1,hp2,hp3,fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
  float tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;

  for (elem=0; elem<actual_size; elem++) {

    for (k=0;k<NGLLZ;k++) {
      for (j=0;j<NGLLY;j++) {
        for (i=0;i<NGLLX;i++) {

          tempx1l = 0.f;
          tempx2l = 0.f;
          tempx3l = 0.f;

          tempy1l = 0.f;
          tempy2l = 0.f;
          tempy3l = 0.f;

          tempz1l = 0.f;
          tempz2l = 0.f;
          tempz3l = 0.f;

          for (l=0;l<NGLLX;l++) {
            hp1 = hprime_xx[l][i];
            tempx1l = tempx1l + dummy_loc[elem][k][j][l][X]*hp1;
            tempy1l = tempy1l + dummy_loc[elem][k][j][l][Y]*hp1;
            tempz1l = tempz1l + dummy_loc[elem][k][j][l][Z]*hp1;

            hp2 = hprime_xx[l][j];
            tempx2l = tempx2l + dummy_loc[elem][k][l][i][X]*hp2;
            tempy2l = tempy2l + dummy_loc[elem][k][l][i][Y]*hp2;
            tempz2l = tempz2l + dummy_loc[elem][k][l][i][Z]*hp2;

            hp3 = hprime_xx[l][k];
            tempx3l = tempx3l + dummy_loc[elem][l][j][i][X]*hp3;
            tempy3l = tempy3l + dummy_loc[elem][l][j][i][Y]*hp3;
            tempz3l = tempz3l + dummy_loc[elem][l][j][i][Z]*hp3;
          }

// compute derivatives of ux, uy and uz with respect to x, y and z
          xixl = jacobian_matrix[elem][k][j][i][XIX];
          xiyl = jacobian_matrix[elem][k][j][i][XIY];
          xizl = jacobian_matrix[elem][k][j][i][XIZ];
          etaxl = jacobian_matrix[elem][k][j][i][ETAX];
          etayl = jacobian_matrix[elem][k][j][i][ETAY];
          etazl = jacobian_matrix[elem][k][j][i][ETAZ];
          gammaxl = jacobian_matrix[elem][k][j][i][GAMMAX];
          gammayl = jacobian_matrix[elem][k][j][i][GAMMAY];
          gammazl = jacobian_matrix[elem][k][j][i][GAMMAZ];
          jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)-xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl));

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl;
          duxdxl_plus_duzdzl = duxdxl + duzdzl;
          duydyl_plus_duzdzl = duydyl + duzdzl;
          duxdyl_plus_duydxl = duxdyl + duydxl;
          duzdxl_plus_duxdzl = duzdxl + duxdzl;
          duzdyl_plus_duydzl = duzdyl + duydzl;

// compute isotropic elements
          kappal = kappa_and_mu[elem][k][j][i][KAPPA];
          mul = kappa_and_mu[elem][k][j][i][MU];

          lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
          lambdal = lambdalplus2mul - 2.f*mul;

// compute stress sigma
          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

          sigma_xy = mul*duxdyl_plus_duydxl;
          sigma_xz = mul*duzdxl_plus_duxdzl;
          sigma_yz = mul*duzdyl_plus_duydzl;

// form dot product with test vector
      tempx1[k][j][i] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl);
      tempy1[k][j][i] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl);
      tempz1[k][j][i] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

      tempx2[k][j][i] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl);
      tempy2[k][j][i] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl);
      tempz2[k][j][i] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

      tempx3[k][j][i] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl);
      tempy3[k][j][i] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl);
      tempz3[k][j][i] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);

        }
      }
    }

    for (k=0;k<NGLLZ;k++) {
      for (j=0;j<NGLLY;j++) {
        for (i=0;i<NGLLX;i++) {

          tempx1l = 0.f;
          tempy1l = 0.f;
          tempz1l = 0.f;

          tempx2l = 0.f;
          tempy2l = 0.f;
          tempz2l = 0.f;

          tempx3l = 0.f;
          tempy3l = 0.f;
          tempz3l = 0.f;

          for (l=0;l<NGLLX;l++) {
            fac1 = hprimewgll_xx[i][l];
            tempx1l = tempx1l + tempx1[k][j][l]*fac1;
            tempy1l = tempy1l + tempy1[k][j][l]*fac1;
            tempz1l = tempz1l + tempz1[k][j][l]*fac1;

            fac2 = hprimewgll_xx[j][l];
            tempx2l = tempx2l + tempx2[k][l][i]*fac2;
            tempy2l = tempy2l + tempy2[k][l][i]*fac2;
            tempz2l = tempz2l + tempz2[k][l][i]*fac2;

            fac3 = hprimewgll_xx[k][l];
            tempx3l = tempx3l + tempx3[l][j][i]*fac3;
            tempy3l = tempy3l + tempy3[l][j][i]*fac3;
            tempz3l = tempz3l + tempz3[l][j][i]*fac3;
          }

          fac1 = wgllwgll_all[k][j][YZ];
          fac2 = wgllwgll_all[k][i][XZ];
          fac3 = wgllwgll_all[j][i][XY];

          sum_terms[elem][k][j][i][X] = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
          sum_terms[elem][k][j][i][Y] = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
          sum_terms[elem][k][j][i][Z] = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

        }
      }
    }

  }

}
#endif // of USE_DEVILLE_INLINED_PRODUCTS

//#pragma omp task input(actual_size,  displ ) inout(max) 
#pragma omp task in(displ[0;actual_size][0;NDIM]) out (*max)
void compute_max(int actual_size,
                 float displ[actual_size][NDIM],
                 float *max)
{
   float current_value;
   int i,j,k,iglob;
   float local_max = -1.f;
//printf("compute_max actual_size %d \n", actual_size);

	for (iglob=0; iglob < actual_size; iglob++){
		current_value = sqrtf(displ[iglob][X]*displ[iglob][X] + displ[iglob][Y]*displ[iglob][Y] + displ[iglob][Z]*displ[iglob][Z]);
                if(current_value > local_max) { local_max = current_value; }
      }
if (local_max> *max){ 
//#pragma css mutex lock (max)
      if (local_max > *max) { *max = local_max; }
//#pragma css mutex unlock (max)
}

}

#pragma omp task out (*seismogram)
void record_seis (int actual_size, float *seismogram, float displ[actual_size][NDIM], int index){
//printf("record_seis actual_size = %d\n", actual_size);
   *seismogram = displ[index][Z];
}

////////////////////////////////////////////////////////////////////////
//                               END   TASKS                          //
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

// global displacement, velocity and acceleration vectors
  static float displ[NGLOB][NDIM];
  static float veloc[NGLOB][NDIM];
  static float accel[NGLOB][NDIM];

// global diagonal mass matrix
  static float rmass_inverse[NGLOB];

// record a seismogram to check that the simulation went well
  static float seismogram[NSTEP];

// arrays with mesh parameters per slice
  static int   ibool[NSPEC][NGLLZ][NGLLY][NGLLX];

  static float jacobian_matrix[NSPEC][NGLLZ][NGLLY][NGLLX][NDIM*NDIM];

  static float kappa_and_mu[NSPEC][NGLLZ][NGLLY][NGLLX][2];

  static float dummy_loc[BS_NSPEC][NGLLZ][NGLLY][NGLLX][NDIM] __attribute__((aligned(16)));

// array with derivatives of Lagrange polynomials and precalculated products
  static float hprime_xx[NGLLX][NGLLX];
  static float hprimewgll_xx[NGLLX][NGLLX];
  static float wgllwgll_all[NGLLY][NGLLX][NDIM];
#ifdef USE_DEVILLE_INLINED_PRODUCTS
  static float hprime_xxT[NGLLX][NGLLX];
  static float hprimewgll_xxT[NGLLX][NGLLX];
  static float all_matrices[NDIM+4][NGLLX][NGLLX];
#endif

  static float sum_terms[BS_NSPEC][NGLLZ][NGLLY][NGLLX][NDIM];

// time step
  int it;

  int ispec,i,j,k,iglob_source;
  int inode, node_size;
  int ii;

// blocking controls
  int actual_size;
  int bs_nspec;
  int bs_nspec_th;

  float Usolidnorm,time,memory_size;

  int i_map_ibool;

// to read external files
  FILE *IIN;

  int connectivity[NSPEC/BS_NSPEC+1][NGLOB/BS_NGLOB+1];
  printf("\nNSPEC = %d\n",NSPEC);
  printf("NGLOB = %d\n\n",NGLOB);
  printf("NSTEP = %d\n",NSTEP);
  printf("deltat = %f\n\n",deltat);

// threads_used = atoi (getenv("CSS_NUM_CPUS"));

#ifdef PROP_THREADS
  int threads_used = 1;
  bs_nspec_th = NSPEC / threads_used;
#else
 bs_nspec_th = BS_NSPEC;
#endif
 
// make sure that we can use the Deville et al. (2002) inlined products
#ifdef USE_DEVILLE_INLINED_PRODUCTS
  if(NGLLX != 5 || NGLLY != 5 || NGLLZ != 5) {
         fprintf(stderr,"we must have NGLLX = NGLLY = NGLLZ = 5 to be able to use the Deville et al. (2002) inlined products, exiting...\n");
         exit(1);
       }
#endif

// estimate total memory size (the size of a real number is 4 bytes)
// we perform the calculation in single precision rather than integer
// to avoid integer overflow in the case of very large meshes
 memory_size = 4.f * ((3.f*NDIM + 1.f) * NGLOB + 13.f * (float)(NGLLX*NGLLY*NGLLZ)*(float)(NSPEC));
 printf("approximate total memory size used = %f Mb\n\n",memory_size/1024.f/1024.f);


// read the mesh from external file
  printf("reading file proc000000_reg1_database.dat\n");
    if((IIN=fopen("proc000000_reg1_database.dat","r"))==NULL) {
          fprintf(stderr,"Cannot open file ./proc000000_reg1_database.dat, exiting...\n");
          exit(1);
        }

  for (ispec=0;ispec<NSPEC;ispec++) {
    for (k=0;k<NGLLZ;k++) {
      for (j=0;j<NGLLY;j++) {
        for (i=0;i<NGLLX;i++) {
// read real numbers here
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][XIX]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][XIY]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][XIZ]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][ETAX]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][ETAY]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][ETAZ]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][GAMMAX]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][GAMMAY]);
          fscanf(IIN, "%e\n", &jacobian_matrix[ispec][k][j][i][GAMMAZ]);
          fscanf(IIN, "%e\n", &kappa_and_mu[ispec][k][j][i][KAPPA]);
          fscanf(IIN, "%e\n", &kappa_and_mu[ispec][k][j][i][MU]);

// read an integer here
          fscanf(IIN, "%d\n", &ibool[ispec][k][j][i]);
// subtract one because indices start at zero in C but this array was created by a Fortran
// program and therefore starts at one in the file stored on the disk
          ibool[ispec][k][j][i]--;
        }
      }
    }
  }
  for (i=0;i<NGLOB;i++) {
    fscanf(IIN, "%e\n", &rmass_inverse[i]);
// the real exactly diagonal mass matrix is read (not its inverse)
// therefore invert it here once and for all
    rmass_inverse[i] = 1.f / rmass_inverse[i];
  }
  fclose(IIN);

// read the derivation matrices from external file
#ifndef USE_PATHS_ROSA
 printf("reading file ./matrices.dat\n");
   if((IIN=fopen("./matrices.dat","r"))==NULL) {
#else
   printf("reading file ./matrices.dat\n");
   if((IIN=fopen("./matrices.dat","r"))==NULL) {

/* printf("reading file ./DB/matrices.dat\n");
   if((IIN=fopen("./DB/matrices.dat","r"))==NULL) {*/
#endif
         fprintf(stderr,"Cannot open file ../multi_GPU_MPI/DATABASES_FOR_SOLVER/matrices.dat, exiting...\n");
         exit(1);
       }

 for (j=0;j<NGLLY;j++) {
   for (i=0;i<NGLLX;i++) {
     fscanf(IIN, "%e\n", &hprime_xx[j][i]);
     fscanf(IIN, "%e\n", &hprimewgll_xx[j][i]);
     fscanf(IIN, "%e\n", &wgllwgll_all[j][i][YZ]);
     fscanf(IIN, "%e\n", &wgllwgll_all[j][i][XZ]);
     fscanf(IIN, "%e\n", &wgllwgll_all[j][i][XY]);

#ifdef USE_DEVILLE_INLINED_PRODUCTS
// DK DK also store the transpose matrices
     hprime_xxT[i][j] = hprime_xx[j][i];
     hprimewgll_xxT[i][j] = hprimewgll_xx[j][i];
#endif

   }
 }

// DK DK store all the matrices in one array in order to reduce the number of arguments sent to the tasks
#ifdef USE_DEVILLE_INLINED_PRODUCTS
 for (j=0;j<NGLLY;j++) {
   for (i=0;i<NGLLX;i++) {
     all_matrices[FLAG_hprime_xx][j][i] = hprime_xx[j][i];
     all_matrices[FLAG_hprime_xxT][j][i] = hprime_xxT[j][i];

     all_matrices[FLAG_hprimewgll_xx][j][i] = hprimewgll_xx[j][i];
     all_matrices[FLAG_hprimewgll_xxT][j][i] = hprimewgll_xxT[j][i];

     all_matrices[FLAG_wgllwgll_YZ][j][i] = wgllwgll_all[j][i][YZ];
     all_matrices[FLAG_wgllwgll_XZ][j][i] = wgllwgll_all[j][i][XZ];
     all_matrices[FLAG_wgllwgll_XY][j][i] = wgllwgll_all[j][i][XY];
   }
 }
#endif

 fclose(IIN);



zones (connectivity, ibool);

//#pragma css start
// clear initial vectors before starting the time loop
// it is not really crucial to parallelize this task because it is done only once
// before entering the time loop therefore it could remain serial
   for (i=0;i<NGLOB;i+=BS_NGLOB) {
      actual_size =  ((NGLOB-i)>=BS_NGLOB ? BS_NGLOB : (NGLOB-i) );
      clear (actual_size, (void*)&displ[i][0], (void*)&veloc[i][0], (void*)&accel[i][0]);
   }

#pragma omp taskwait
  printf("\nstarting the time loop\n\n");
  t_start = usecs();

// start of the time loop (which must remain serial obviously)
  for (it=1;it<=NSTEP;it++) {

// compute maximum of norm of displacement from time to time and display it
// in order to monitor the simulation
// this can remain serial because it is done only every NTSTEP_BETWEEN_OUTPUT_INFO time steps
   if((it % NTSTEP_BETWEEN_OUTPUT_INFO) == 0 || it == 5 || it == NSTEP) {

      Usolidnorm = -1.f;
      for (i=0;i<NGLOB;i+=BS_NGLOB) {
        actual_size =  ((NGLOB-i)>=BS_NGLOB ? BS_NGLOB : (NGLOB-i) );
        compute_max(actual_size,  (void*)&displ[i][0],  &Usolidnorm);
       }

//#pragma css wait on (&Usolidnorm)
#pragma omp taskwait
      printf("Time step # %d out of %d\n",it,NSTEP);
      printf("Max norm displacement vector U in the solid (m) = %.8g\n\n",Usolidnorm);

// check stability of the code, exit if unstable
    if(Usolidnorm > STABILITY_THRESHOLD) {
       fprintf(stderr,"code became unstable and blew up\n");
       exit(1);
      }
    }

//#pragma css barrier
// big loop over all the global points (not elements) in the mesh to update
// the displacement and velocity vectors and clear the acceleration vector
    for (i=0;i<NGLOB;i+=BS_NGLOB) {
      actual_size =  ((NGLOB-i)>=BS_NGLOB ? BS_NGLOB : (NGLOB-i) );
      update_disp_vel (actual_size, (void*)&displ[i][0], (void*)&veloc[i][0], (void*)&accel[i][0]);
    }

//#pragma css barrier
// BArrier afegit perque dona calcul incorrecte de usolidnorm

// big loop over all the elements in the mesh to localize data once and for all
// from the global vectors to the local mesh
// using indirect addressing (contained in array ibool)

// big loop over all the elements in the mesh to compute the elemental contribution
// to the acceleration vector of each element of the finite-element mesh

// need to ensure ALL displ has been produced.
// Limitation of current dependence detection mechanism.
// this is waiting for displ[] to be entirely filled because we use it in the gather below.
//#pragma css barrier
//    for (ispec=0;ispec<NSPEC;ispec+=BS_NSPEC) {
    bs_nspec = bs_nspec_th;
    for (ispec=0;ispec<NSPEC;ispec+=bs_nspec) {

//      actual_size =  ((NSPEC-ispec)>=BS_NSPEC ? BS_NSPEC : (NSPEC-ispec) );
      actual_size =  ((NSPEC-ispec)>=bs_nspec ? bs_nspec : (NSPEC-ispec) );
// DK DK to Rosa: I put the (void*) back in front of &ibool otherwise I get a warning
// DK DK from both GNU gcc and Intel icc again

       for (inode=0; inode < NGLOB; inode+=BS_NGLOB)  {
         node_size =  ((NGLOB-inode)>=BS_NGLOB ? BS_NGLOB : (NGLOB-inode) );
	 if (connectivity[ispec/BS_NSPEC][inode/BS_NGLOB]!=0){
      		gather (actual_size, node_size, (void*)&displ[inode][0], (void*)&ibool[ispec][0][0][0], dummy_loc, inode/BS_NGLOB );
	}	

       }
//#pragma css barrier
      process_element(actual_size, dummy_loc,
#ifdef USE_DEVILLE_INLINED_PRODUCTS
         all_matrices,
#else
         hprime_xx,
         hprimewgll_xx,
         wgllwgll_all,
#endif
         (void*)&jacobian_matrix[ispec][0][0][0][0],
         (void*)&kappa_and_mu[ispec][0][0][0][0],
         sum_terms );

// sum contributions from each element to the global mesh using indirect addressing
//#pragma css barrier
       for (inode=0; inode < NGLOB; inode+=BS_NGLOB)  {
         node_size =  ((NGLOB-inode)>=BS_NGLOB ? BS_NGLOB : (NGLOB-inode) );
	 if (connectivity[ispec/BS_NSPEC][inode/BS_NGLOB]!=0){
//		printf("gather node = %d zona = %d \n", ispec/BS_NSPEC, inode/BS_NGLOB);
      //		gather (actual_size, node_size, &displ[inode][0], &ibool[ispec][0][0][0], dummy_loc, inode/BS_NGLOB );
      //scatter (actual_size, &ibool[ispec][0][0][0], sum_terms, accel);
      		scatter (actual_size, node_size,  (void*)&ibool[ispec][0][0][0], sum_terms, (void*)&accel[inode][0], inode/BS_NGLOB);
	}	
      }


/*      if ((bs_nspec==bs_nspec_th)&&((NSPEC-ispec) <= (NSPEC/10))) {
	bs_nspec = bs_nspec >> 1;
      }*/
    }   // end of main loop on all the elements

//#pragma css barrier

// add the earthquake source at a given grid point
// this is negligible and is intrinsically serial because it is done by only
// one grid point out of several millions typically
// we subtract one to the element number of the source because arrays start at 0 in C
// compute current time
  time = (it-1)*deltat;
  iglob_source = ibool[NSPEC_SOURCE-1][1][1][1];
// we divide the amplitude of the source by rmass_inverse[iglob_source] here because
// we have merged the calculation of acceleration and velocity below in a single task
// and therefore the value of accel[] at that point will be
// multiplied by rmass_inverse[i] in that merged task
//#pragma css barrier
 // accel[iglob_source][Z] += 1.e4f * (1.f - 2.f*a*(time-t0)*(time-t0)) * expf(-a*(time-t0)*(time-t0)) / (rho * rmass_inverse[iglob_source]);
//This line has been added to update_acc_vel

// big loop over all the global points (not elements) in the mesh to update
// the acceleration and velocity vectors
    for (i=0;i<NGLOB;i+=BS_NGLOB) {
      actual_size =  ((NGLOB-i)>=BS_NGLOB ? BS_NGLOB : (NGLOB-i) );
      if ((iglob_source >= i) && (iglob_source < i + actual_size))
         update_acc_vel (actual_size, (void*)&accel[i][0], (void*)&veloc[i][0], (void*)&rmass_inverse[i], iglob_source%BS_NGLOB, time);
	
      else update_acc_vel (actual_size, (void*)&accel[i][0], (void*)&veloc[i][0], (void*)&rmass_inverse[i], -1, 0.0); 
    }

//ULTIM
//#pragma css barrier
// record a seismogram to check that the simulation went well
// we subtract one to the element number of the receiver because arrays start at 0 in C
//   seismogram[it-1] = displ[ibool[NSPEC_STATION-1][1][1][1]][Z];
      int chunck = (ibool[NSPEC_STATION-1][1][1][1]/BS_NGLOB)*BS_NGLOB;
      actual_size =  ((NGLOB-chunck)>=BS_NGLOB ? BS_NGLOB : (NGLOB-chunck) );
//printf ("actual size %d\n", actual_size);
//printf("ibool= %d, index = %d, sector=%d\n", ibool[NSPEC_STATION-1][1][1][1], ibool[NSPEC_STATION-1][1][1][1]- (ibool[NSPEC_STATION-1][1][1][1]/BS_NGLOB)*BS_NGLOB, (ibool[NSPEC_STATION-1][1][1][1]/BS_NGLOB)*BS_NGLOB);
     record_seis (actual_size, &seismogram[it-1], &displ[chunck][0], ibool[NSPEC_STATION-1][1][1][1] - chunck);

  } // end of the serial time loop

#pragma omp taskwait
  t_end = usecs();

//#pragma css finish
  printf("elapsed time: %f seconds\n",(float)(t_end-t_start)/1000000.f);

// save the seismogram at the end of the run
// vdimic: removed writing to file, so it doesn't get traced
// char filename[50];
// sprintf(filename, "seismogram_SMP_StarSs_%d.txt\0",getpid());
// if((IIN = fopen(filename,"w")) == NULL) {
//         fprintf(stderr,"Cannot open, exiting...\n");
//         exit(1);
//       }
// for (it=0;it<NSTEP;it++)
// {  fprintf(IIN,"%e %e\n",it*deltat,seismogram[it]);
// }
// fclose(IIN);

   printf("\nEnd of the program\n\n");

}

//
// function to measure time
//

long usecs (){
  struct timeval t;

  gettimeofday(&t,NULL);
  return t.tv_sec*1000000+t.tv_usec;
}

