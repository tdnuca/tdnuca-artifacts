#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <sys/time.h>
//#include <time.h>

//using enums to get real C constants
//because defines cannot be used for OmpSs pragmas
enum { NB=24 };
enum { BS=128 };
//enum { NB=16 };
//enum { BS=64 };
//enum { NB=8 };
//enum { BS=8 };

#define FLOATS_PER_LINE 15

#define FALSE (0)
#define TRUE (1)


#ifdef SPU_CODE
#include <spu_intrinsics.h>
#endif

typedef float (*p_block_t)[BS];

p_block_t A[NB][NB];
p_block_t origA[NB][NB];
p_block_t L[NB][NB];
p_block_t U[NB][NB];


p_block_t allocate_clean_block()
{
  int i,j;
  p_block_t p;

  p=(p_block_t)malloc(BS*BS*sizeof(float));
  if (p!=NULL){
     for (i = 0; i < BS; i++) 
        for (j = 0; j < BS; j++) 
           p[i][j]=0.0; 
  }
  else {
     printf ("OUT OF MEMORY!!!!!!!!!!!!!!!\n");
     exit (-1);
  }
  return (p);
}


void genmat ()
{
   int  i, j, ii, jj;
   int null_entry;

   long int init_val;

   init_val = 1325;

/* structure */
   for (ii=0; ii<NB; ii++)
      for (jj=0; jj<NB; jj++){
         null_entry=FALSE;
         if ((ii<jj) && (ii%3 !=0)) null_entry =TRUE;
         if ((ii>jj) && (jj%3 !=0)) null_entry =TRUE;
	 if (ii%2==1) null_entry=TRUE;
	 if (jj%2==1) null_entry=TRUE;
	 if (ii==jj) null_entry=FALSE;
	 if (ii==jj-1) null_entry=FALSE;
         if (ii-1 == jj) null_entry=FALSE; 
         if (null_entry==FALSE){
            A[ii][jj] = allocate_clean_block();
	    if (A[ii][jj]==NULL) {
		printf("Out of memory\n");
		exit(1);
	    }
         }
         else A[ii][jj] = NULL;
      }

   /* Initialization */
   for (ii = 0; ii < NB; ii++)
      for (jj = 0; jj < NB; jj++)
      {
         p_block_t p;
         p = A[ii][jj];
         if (p!=NULL)
         for (i = 0; i < BS; i++)
            for (j = 0; j < BS; j++) {
               init_val = (3125 * init_val) % 65536;
               p[i][j] = 0.0001;
               if (ii == jj){
                  if (i==j) p[i][j] = -20000 ;
                     if (((i-1)==j) || (i==(j-1))) p[i][j] = 10000 ;
               }
         }
   }
}


#pragma omp task in (ref_block[0;BS][0;BS], to_comp[0;BS][0;BS]) out (*mse)
void are_blocks_equal (float ref_block[BS][BS], float to_comp[BS][BS], float *mse)
{
   int i,j;
   float diff;

   *mse = 0.0;
   for (i = 0; i < BS; i++) 
     for (j = 0; j < BS; j++) {
        diff = ref_block[i][j]-to_comp[i][j];
        *mse += diff*diff;
     }
}

void compare_mat (p_block_t X[NB][NB], p_block_t Y[NB][NB], struct timeval *stop)
{
   int ii, jj;
   float sq_error[NB][NB];
   p_block_t Zero_block;
   int some_difference = FALSE;

   Zero_block = allocate_clean_block();

   for (ii = 0; ii < NB; ii++) 
     for (jj = 0; jj < NB; jj++) {
       if (X[ii][jj] == NULL) 
          if (Y[ii][jj] == NULL)
              sq_error[ii][jj] = 0.0f;
          else
              are_blocks_equal(Zero_block, Y[ii][jj],&sq_error[ii][jj]);
       else
          are_blocks_equal(X[ii][jj], Y[ii][jj],&sq_error[ii][jj]);
     }

/* Alternative to put wait on */
#pragma omp taskwait 
gettimeofday(stop,NULL);
//   printf ("\nComparison of matrices at %x and %x\n",X,Y);
   for (ii = 0; ii < NB; ii++) 
     for (jj = 0; jj < NB; jj++) 
         if (sq_error[ii][jj] >0.0000001L) {
 //            printf ("block [%d, %d]: detected mse = %.20lf\n",ii,jj,sq_error[ii][jj]);
              some_difference =TRUE;
             // exit(1);
         }
     if ( some_difference == FALSE) printf ("matrices are identical\n");

}


#pragma omp task inout(diag[0;BS][0;BS])
void lu0(float diag[BS][BS])
{
   int i, j, k;

   for (k=0; k<BS; k++)
      for (i=k+1; i<BS; i++) {
         diag[i][k] = diag[i][k] / diag[k][k];
         for (j=k+1; j<BS; j++)
             diag[i][j] -= diag[i][k] * diag[k][j];
      }

}

#pragma omp task in(diag[0;BS][0;BS]) inout(row[0;BS][0;BS])
void bdiv(float diag[BS][BS], float row[BS][BS])
{
   int i, j, k;

   for (i=0; i<BS; i++)
      for (k=0; k<BS; k++) {
         row[i][k] = row[i][k] / diag[k][k];
         for (j=k+1; j<BS; j++)
            row[i][j] -= row[i][k]*diag[k][j];
      }

}


#pragma omp task in(row[0;BS][0;BS],col[0;BS][0;BS]) inout(inner[0;BS][0;BS])
void bmod(float row[BS][BS], float col[BS][BS], float inner[BS][BS])
{
  int i, j, k;


  for (i=0; i<BS; i++)
     for (j=0; j<BS; j++)
        for (k=0; k<BS; k++) 
           inner[i][j] -= row[i][k]*col[k][j];

}

#pragma omp task in(a[0;BS][0;BS],b[0;BS][0;BS]) inout(c[0;BS][0;BS])
void block_mpy_add(float a[BS][BS], float b[BS][BS], float c[BS][BS])
{
  int i, j, k;

  for (i=0; i<BS; i++)
     for (j=0; j<BS; j++)
        for (k=0; k<BS; k++) 
           c[i][j] += a[i][k]*b[k][j];
}


#pragma omp task in(diag[0;BS][0;BS]) inout(col[0;BS][0;BS])
void fwd(float diag[BS][BS], float col[BS][BS])
{
  int i, j, k;

  for (j=0; j<BS; j++)
     for (k=0; k<BS; k++) 
        for (i=k+1; i<BS; i++)
           col[i][j] -= diag[i][k]*col[k][j];

}

#pragma omp task in(A[0;BS][0;BS]) out(L[0;BS][0;BS], U[0;BS][0;BS])
void split_block (float A[BS][BS], float L[BS][BS], float U[BS][BS]) 
{
  int i, j, k;

  for (i=0; i<BS; i++)
     for (j=0; j<BS; j++) {
        if (i==j)     { L[i][j]=1.0;      U[i][j]=A[i][j]; }
        else if (i>j) { L[i][j]=A[i][j];  U[i][j]=0.0;}
       else           { L[i][j]=0.0;      U[i][j]=A[i][j]; }
     }
}


#pragma omp task in(Src[0;BS][0;BS]) out(Dst[0;BS][0;BS])
void copy_block (float Src[BS][BS], float Dst[BS][BS])
{
  int i, j;

  for (i=0; i<BS; i++)
     for (j=0; j<BS; j++) 
        Dst[i][j]=Src[i][j];
}
#pragma omp task out(Dst[0;BS][0;BS])
void clean_block (float Dst[BS][BS] )
{
  int i, j;

  for (i=0; i<BS; i++)
     for (j=0; j<BS; j++) 
        Dst[i][j]=0.0;
}

void LU(p_block_t A[NB][NB])
{
   int ii, jj, kk;

   for (kk=0; kk<NB; kk++) {
      lu0( A[kk][kk]);
      for (jj=kk+1; jj<NB; jj++)
         if (A[kk][jj] != NULL) fwd(A[kk][kk], A[kk][jj]);

      for (ii=kk+1; ii<NB; ii++) 
         if (A[ii][kk] != NULL)
            bdiv (A[kk][kk], A[ii][kk]);

      for (ii=kk+1; ii<NB; ii++) {
         if (A[ii][kk] != NULL) {
            for (jj=kk+1; jj<NB; jj++) {
               if (A[kk][jj] != NULL) {
                  if (A[ii][jj]==NULL) A[ii][jj]=allocate_clean_block();
                  bmod(A[ii][kk], A[kk][jj], A[ii][jj]);
               }
            }
         }
      }
   }
}

void split_mat (p_block_t LU[NB][NB],p_block_t L[NB][NB],p_block_t U[NB][NB])
{
  int ii, jj;
  p_block_t block;

  for (ii=0; ii<NB; ii++) 
     for (jj=0; jj<NB; jj++){
        if (ii==jj) {                               /* split diagonal block */
           L[ii][ii] = allocate_clean_block(); 
           U[ii][ii] = allocate_clean_block(); 
           split_block (LU[ii][ii],L[ii][ii],U[ii][ii]);
        } else {                                    /* copy non diagonal block to ... */
            if (LU[ii][jj] != NULL) {
              block = allocate_clean_block(); 
              copy_block(LU[ii][jj],block);
           } else block = NULL;
           if (ii>jj) {                             /*         ...either L ... */
              L[ii][jj]=block;
              U[ii][jj]=NULL;
           } else {                                /*         ... or U */
              L[ii][jj]=NULL;
              U[ii][jj]=block;
          }
        }
      }
}

void copy_mat (p_block_t Src[NB][NB], p_block_t Dst[NB][NB])
{
  int ii, jj;
  p_block_t block;

  for (ii=0; ii<NB; ii++)
     for (jj=0; jj<NB; jj++)
        if (Src[ii][jj] != NULL) {
           block = allocate_clean_block();
           copy_block(Src[ii][jj],block);
           Dst[ii][jj] = block;
        } else
           Dst[ii][jj]=NULL;
}

void clean_mat (p_block_t Src[NB][NB])
{
  int ii, jj;

  for (ii=0; ii<NB; ii++)
     for (jj=0; jj<NB; jj++)
        if (Src[ii][jj] != NULL) {
/* this will work in sequential, but would require waiting for all uses of Src[ii][jj] in parallel
           free (Src[ii][jj]);
           Src[ii][jj]=NULL;
*/
        clean_block(Src[ii][jj]);
        }
}

/* C = A*B */
void sparse_matmult (p_block_t A[NB][NB], p_block_t B[NB][NB], p_block_t C[NB][NB])
{
  int ii, jj, kk;

  for (ii=0; ii<NB; ii++) 
     for (jj=0; jj<NB; jj++)
        for (kk=0; kk<NB; kk++)
           if ((A[ii][kk] != NULL) && (B[kk][jj] !=NULL )) {
              if (C[ii][jj] == NULL) C[ii][jj] = allocate_clean_block();
              block_mpy_add (A[ii][kk], B[kk][jj], C[ii][jj]);
           }
}

int main(int argc, char* argv[])
{

  struct timeval start;
  struct timeval stop;
  unsigned long elapsed;
 
   


   genmat();
   gettimeofday(&start,NULL);

   copy_mat (A, origA); 

   LU (A);

   split_mat (A, L, U);

   clean_mat (A);

   sparse_matmult (L, U, A); 
   compare_mat (origA, A, &stop);

#pragma omp taskwait

  gettimeofday(&stop,NULL);
  elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
  elapsed += stop.tv_usec - start.tv_usec;
  //printf ("Time %lu microsecs\n", elapsed);
  printf ("par_sec_time_us:%lu\n", elapsed);

}

