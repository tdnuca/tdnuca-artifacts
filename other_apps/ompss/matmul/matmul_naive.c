/*
* Copyright (c) 2008, BSC (Barcelon Supercomputing Center)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY BSC ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

//const int BSIZE = 64; original
const int BSIZE = 128;

//-------------------------------------------------------------------------------------------------------------------------

// NB: could also use reduction(+:[BSIZE][BSIZE]C)
#pragma omp task in(A[0;BSIZE][0;BSIZE], B[0;BSIZE][0;BSIZE]) commutative(C[0;BSIZE][0;BSIZE])
void matmul (float A[BSIZE][BSIZE], float B[BSIZE][BSIZE], float C[BSIZE][BSIZE])
{
  int i, j, k;

  for (i = 0; i < BSIZE; i++)
  {
    for (j = 0; j < BSIZE; j++)
	{
      float c = 0.0;
      for (k=0; k < BSIZE; k++)
        c += A[i][k]*B[k][j];
      C[i][j] += c;
    }
  }
}

//---------------------------------------------------------------------------------------------------------------------------

void compute(struct timeval *start, struct timeval *stop, int DIM, float *A[DIM][DIM], float *B[DIM][DIM], float *C[DIM][DIM])
{
 int i, j, k;

  gettimeofday(start,NULL);

  for (i = 0; i < DIM; i++)
    for (j = 0; j < DIM; j++)
      for (k = 0; k < DIM; k++){
        matmul ( (float (*)[BSIZE])A[i][k], (float (*)[BSIZE])B[k][j], (float (*)[BSIZE])C[i][j]); //arico: cast the values to match the matmul function arguments
}
#pragma omp taskwait
  gettimeofday(stop,NULL);

}

float **A;
float **B;
float **C;

void initialize (int argc, char **argv, int * N_p, int * DIM_p);

int main (int argc, char **argv)
{
  // local vars
  int i, j, k;
  int N,DIM;
  struct timeval start;
  struct timeval stop;
  unsigned long elapsed;


  // application inicializations
  initialize (argc, argv, &N, &DIM);

  // compute with CellSs
  compute(&start, &stop, DIM, (void *)A, (void *)B, (void *)C);

  elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
  elapsed += stop.tv_usec - start.tv_usec;
  printf("Matrix dimension: %d\n",N);
  printf ("Time %lu microsecs\n", elapsed);
  printf ("Perf %lu MFlops\n", (long unsigned)(((double)N*N*N*2)/((double)elapsed)));

  printf("par_sec_time_us:%lu\n",elapsed);

  return 0;
}

static void convert_to_blocks(int DIM, int N, float *Alin, float *A[DIM][DIM])
{
  int i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      A[i/BSIZE][j/BSIZE][(i%BSIZE)*BSIZE+j%BSIZE] = Alin[j*N+i];
    }
  }

}

void fillRandom(float *Alin, int NN)
{
    int i;
    for (i = 0; i < NN; i++)
    {
        Alin[i]=((float)rand())/((float)RAND_MAX);
    }
}

void initialize (int argc, char **argv, int * N_p, int * DIM_p)
{
  int DIM;
  int i;

  if (argc==2)
  {
    DIM=atoi(argv[1]);
  }
  else
  {
    printf("usage: %s DIM\n",argv[0]);
    exit(0);
  }

  // matrix init
  long int N=BSIZE*DIM;
  long int NN=N*N;

  *N_p=N;
  *DIM_p=DIM;

  // linear matrix
  float *Alin = (float *) malloc(NN * sizeof(float));
  float *Blin = (float *) malloc(NN * sizeof(float));
  float *Clin = (float *) malloc(NN * sizeof(float));

  // fill the matrix with random values
  fillRandom(Alin, NN);
  fillRandom(Blin, NN);
  memset(Clin, 0, NN);

  A = (float **) malloc(DIM*DIM*sizeof(float *));
  B = (float **) malloc(DIM*DIM*sizeof(float *));
  C = (float **) malloc(DIM*DIM*sizeof(float *));

  for (i = 0; i < DIM*DIM; i++)
  {
     A[i] = (float *) malloc(BSIZE*BSIZE*sizeof(float));
     B[i] = (float *) malloc(BSIZE*BSIZE*sizeof(float));
     C[i] = (float *) malloc(BSIZE*BSIZE*sizeof(float));
  }

  convert_to_blocks(DIM, N, Alin, (void *)A);
  convert_to_blocks(DIM, N, Blin, (void *)B);
  convert_to_blocks(DIM, N, Clin, (void *)C);

  free(Alin);
  free(Blin);
  free(Clin);
}

