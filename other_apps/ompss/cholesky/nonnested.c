#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <errno.h>
//#include <cuda_runtime.h>

// Configuration options
#define CONVERT_TASK 0
#define CONVERT_REQUEST 0
#define POTRF_SMP 1
#define POTRF_NESTED 0
#define TRSM_SMP 1
#define TRSM_NESTED 0
#define SYRK_SMP 1
#define SYRK_NESTED 0
#define GEMM_SMP 1
#define GEMM_NESTED 0
#define MAGMA_BLAS 0
#define USE_AFFINITY 0
#define MKL 0 
//#define DOUBLE_PREC 1
#define ATLAS 1

#define STOP_SCHED 0 
#define DELAY 1000000 
#define CONDITION 0
//( sched_getcpu() != 1 && sched_getcpu() != 2 )
//&& sched_getcpu() != 3 && sched_getcpu() != 4 && sched_getcpu() != 5 && sched_getcpu() != 6 && sched_getcpu() != 7 && sched_getcpu() != 8 )
//&& sched_getcpu() != 9 && sched_getcpu() != 10 && sched_getcpu() != 11 && sched_getcpu() != 12 && sched_getcpu() != 13 && sched_getcpu() != 14 && sched_getcpu() != 15 && sched_getcpu() != 16)
int NUM_NODES;

#if !(POTRF_SMP==1 && TRSM_SMP==1 && SYRK_SMP==1 && GEMM_SMP==1 )
#  include <cuda_runtime.h>
#endif
// Checking constraints
// If converting matrix under request, mark convert functions as tasks
#if CONVERT_REQUEST
 #if !CONVERT_TASK
  #undef CONVERT_TASK
  #define CONVERT_TASK 1
 #endif
#endif

// Include GPU kernel's library: MAGMA or CUBLAS
#if MAGMA_BLAS
#include <magma.h>
#else
#if MKL
#include <mkl.h>
#else
#if ATLAS
#include <clapack.h>
#include <cblas.h>
extern void spotrf_(char *, int *, float *, int *, int *);
#else
#include <cublas.h>
#endif
#endif
#endif 
// Define macro's to make the code cleaner
#if POTRF_NESTED
 #define CALL_POTRF_TILE(x1, x2) smp_cholesky(x1, x2);
#else
 #define CALL_POTRF_TILE(x1, x2) potrf_tile(x1, x2);
#endif

// Define macro's to make the code cleaner
#if TRSM_NESTED
#define CALL_TRSM_TILE(x1, x2, x3) omp_trsm(x1, x2, x3);
#else
#define CALL_TRSM_TILE(x1, x2, x3) trsm_tile(x1, x2, x3);
#endif

#if SYRK_NESTED
#define CALL_SYRK_TILE(x1, x2, x3) omp_syrk(x1, x2, x3);
#else
#define CALL_SYRK_TILE(x1, x2, x3) syrk_tile(x1, x2, x3);
#endif

#if GEMM_NESTED
#define CALL_GEMM_TILE(x1, x2, x3, x4) omp_gemm(x1, x2, x3, x4);
#else
#define CALL_GEMM_TILE(x1, x2, x3, x4) gemm_tile(x1, x2, x3, x4);
#endif

#if CONVERT_REQUEST
 #define TRY_GATHER_BLOCK(x1, x2, x3, x4) \
    if (x4 == NULL) { \
        x4 = malloc_block(x2); \
        gather_block(x1, x2, x3, x4); \
    }

 #define CHECK_BLOCK_NOT_NULL(x1) if (A[i][j] != NULL)
#else
 #define TRY_GATHER_BLOCK(x1, x2, x3, x4) 
 #define CHECK_BLOCK_NOT_NULL(x1) 
#endif

#ifdef DOUBLE_PREC
#define REAL double
#define laset_		dlaset_
#define gemm_		dgemm_
#define trsm_       dtrsm_
#define trmm_       dtrmm_
#define syrk_       dsyrk_
#define lange_		dlange_
#define lacpy_		dlacpy_
#define larnv_		dlarnv_
#define potrf_		dpotrf_
#define gpu_potrf       gpu_d_potrf
#define gpu_blas_gemm	gpu_blas_d_gemm
#define gpu_blas_trsm	gpu_blas_d_trsm
#define gpu_blas_syrk	gpu_blas_d_syrk
#define accepted_error	1.0e-15
#else
#define REAL float
#define laset_		slaset_
#define gemm_		sgemm_
#define trsm_       strsm_
#define trmm_       strmm_
#define syrk_       ssyrk_
#define lange_		slange_
#define lacpy_		slacpy_
#define larnv_		slarnv_
#define potrf_		spotrf_
#define gpu_potrf       gpu_s_potrf
#define gpu_blas_gemm	gpu_blas_s_gemm
#define gpu_blas_trsm	gpu_blas_s_trsm
#define gpu_blas_syrk	gpu_blas_s_syrk
#define accepted_error	1.0e-1
#endif


#if MAGMA_BLAS
#define gpu_d_potrf     magma_dpotrf_gpu
#define gpu_blas_d_gemm magmablas_dgemm
#define gpu_blas_d_trsm magmablas_dtrsm
#define gpu_blas_d_syrk magmablas_dsyrk
#define gpu_s_potrf     magma_spotrf_gpu
#define gpu_blas_s_gemm magmablas_sgemm
#define gpu_blas_s_trsm magmablas_strsm
#define gpu_blas_s_syrk cublasSsyrk // = magmablas_ssyrk
#else
#define gpu_blas_d_gemm cublasDgemm
#define gpu_blas_d_trsm cublasDtrsm
#define gpu_blas_d_syrk cublasDsyrk
#define gpu_blas_s_gemm cublasSgemm
#define gpu_blas_s_trsm cublasStrsm
#define gpu_blas_s_syrk cublasSsyrk
#endif

#if MKL
#define LONG int
#else
#define LONG long
#endif

//#include <papi.h>
void laset_ (const char * UPLO, const int * M, const int * N, const REAL * ALPHA, const REAL * BETA, REAL * A, const int * LDA);

void gemm_ (const char *transa, const char *transb, int *l, int *n, int *m, REAL *alpha,
             const void *a, int *lda, void *b, int *ldb, REAL *beta, void *c, int *ldc);

REAL lange_ (const char *norm, const int *m, const int *n, const REAL *a, const int *lda, REAL *work);

void lacpy_ (const char *norm, const int *m, const int *n, const REAL *a, const int *lda, REAL *b, const int *ldb);



void larnv_ (const int *idist, int *iseed, const int *n, REAL *x);


//void potrf_( const char* uplo, const int* n, REAL* a, const int* lda, LONG* info );


void trsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, REAL *alpha, REAL *a, int *lda, REAL *b, int *ldb);
void trmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, REAL *alpha, REAL *a, int *lda, REAL *b, int *ldb);
void syrk_(char *uplo, char *trans, int *n, int *k, REAL *alpha, REAL *a, int *lda, REAL *beta, REAL *c, int *ldc);
float get_time();
static int check_factorization(int, REAL*, REAL*, int, char , REAL);

enum blas_order_type {
            blas_rowmajor = 101,
            blas_colmajor = 102 };

enum blas_cmach_type {
            blas_base      = 151,
            blas_t         = 152,
            blas_rnd       = 153,
            blas_ieee      = 154,
            blas_emin      = 155,
            blas_emax      = 156,
            blas_eps       = 157,
            blas_prec      = 158,
            blas_underflow = 159,
            blas_overflow  = 160,
            blas_sfmin     = 161};

enum blas_norm_type {
            blas_one_norm       = 171,
            blas_real_one_norm  = 172,
            blas_two_norm       = 173,
            blas_frobenius_norm = 174,
            blas_inf_norm       = 175,
            blas_real_inf_norm  = 176,
            blas_max_norm       = 177,
            blas_real_max_norm  = 178 };

static void
BLAS_error(char *rname, int err, int val, int x) {
  fprintf( stderr, "%s %d %d %d\n", rname, err, val, x );
  abort();
}

static
void
BLAS_ge_norm(enum blas_order_type order, enum blas_norm_type norm,
  int m, int n, const REAL *a, int lda, REAL *res) {
  int i, j; float anorm, v;
  char rname[] = "BLAS_ge_norm";

  if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

  if (norm == blas_frobenius_norm) {
    anorm = 0.0f;
    for (j = n; j; --j) {
      for (i = m; i; --i) {
        v = a[0];
        anorm += v * v;
        a++;
      }
      a += lda - m;
    }
    anorm = sqrt( anorm );
  } else if (norm == blas_inf_norm) {
    anorm = 0.0f;
    for (i = 0; i < m; ++i) {
      v = 0.0f;
      for (j = 0; j < n; ++j) {
        v += abs( a[i + j * lda] );
      }
      if (v > anorm)
        anorm = v;
    }
  } else {
    BLAS_error( rname, -2, norm, 0 );
    return;
  }

  if (res) *res = anorm;
}

static
REAL
BLAS_dpow_di(REAL x, int n) {
  REAL rv = 1.0;

  if (n < 0) {
    n = -n;
    x = 1.0 / x;
  }

  for (; n; n >>= 1, x *= x) {
    if (n & 1)
      rv *= x;
  }

  return rv;
}

static
REAL
BLAS_dfpinfo(enum blas_cmach_type cmach) {
  REAL eps = 1.0, r = 1.0, o = 1.0, b = 2.0;
  int t = 53, l = 1024, m = -1021;
  char rname[] = "BLAS_dfpinfo";

  if ((sizeof eps) == sizeof(float)) {
    t = 24;
    l = 128;
    m = -125;
  } else {
    t = 53;
    l = 1024;
    m = -1021;
  }

  /* for (i = 0; i < t; ++i) eps *= half; */
  eps = BLAS_dpow_di( b, -t );
  /* for (i = 0; i >= m; --i) r *= half; */
  r = BLAS_dpow_di( b, m-1 );

  o -= eps;
  /* for (i = 0; i < l; ++i) o *= b; */
  o = (o * BLAS_dpow_di( b, l-1 )) * b;

  switch (cmach) {
    case blas_eps: return eps;
    case blas_sfmin: return r;
    default:
      BLAS_error( rname, -1, cmach, 0 );
      break;
  }
  return 0.0;
}



int ISEED[4] = {0,0,0,1};
int IINFO;
int intONE=1;

//void gpu_spotf2_var1_( char *, int*, unsigned int*, int *, int * );
void gpu_spotrf_var1_( char *, int*, unsigned int*, int *, int *, int * );

void cholesky(REAL *Alin, REAL** Ah, int ts, int nt);

void add_to_diag_hierarchical (REAL ** matrix, int ts, int nt, float alpha)
{
    int i;

    for (i = 0; i < nt * ts; i++)
        matrix[(i/ts) * nt + (i/ts)][(i%ts) * ts + (i%ts)] += alpha;
}

void add_to_diag(REAL * matrix, int n, REAL alpha)
{
    int i;

    for (i = 0; i < n; i++)
        matrix[ i + i * n ] += alpha;
}

double gtod_ref_time_sec = 0.0;

float get_time()
{
    double t, norm_sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);

    // If this is the first invocation of through dclock(), then initialize the
    // "reference time" global variable to the seconds field of the tv struct.
    if (gtod_ref_time_sec == 0.0)
        gtod_ref_time_sec = (double) tv.tv_sec;

    // Normalize the seconds field of the tv struct so that it is relative to the
    // "reference time" that was recorded during the first invocation of dclock().
    norm_sec = (double) tv.tv_sec - gtod_ref_time_sec;

    // Compute the number of seconds since the reference time.
    t = norm_sec + tv.tv_usec * 1.0e-6;

    return (float) t;
}
/*------------------------------------------------------------------------
 *  *  Robust Check the factorization of the matrix A2
 *   */
static int check_factorization(int N, REAL *A1, REAL *A2, int LDA, char uplo, REAL eps)
{
    REAL Anorm, Rnorm;
    REAL alpha;
    int info_factorization;
    int i,j;
    char NORM = 'I', ALL = 'A', UP = 'U', LO = 'L', TR = 'T', NU = 'N', RI = 'R';

    REAL *Residual = (REAL *)malloc(N*N*sizeof(REAL));
    REAL *L1       = (REAL *)malloc(N*N*sizeof(REAL));
    REAL *L2       = (REAL *)malloc(N*N*sizeof(REAL));
    REAL *work     = (REAL *)malloc(N*sizeof(REAL));

    memset((void*)L1, 0, N*N*sizeof(REAL));
    memset((void*)L2, 0, N*N*sizeof(REAL));

    alpha= 1.0;

    lacpy_(&ALL, &N, &N, A1, &LDA, Residual, &N);

    /* Dealing with L'L or U'U  */
    if (uplo == 'U'){
        lacpy_(&UP, &N, &N, A2, &LDA, L1, &N);
        lacpy_(&UP, &N, &N, A2, &LDA, L2, &N);
        trmm_(&LO, &uplo, &TR, &NU, &N, &N, &alpha, L1, &N, L2, &N);
    }
    else{
        lacpy_(&LO, &N, &N, A2, &LDA, L1, &N);
        lacpy_(&LO, &N, &N, A2, &LDA, L2, &N);
        trmm_(&RI, &LO, &TR, &NU, &N, &N, &alpha, L1, &N, L2, &N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    Rnorm = lange_(&NORM, &N, &N, Residual, &N, work);
    Anorm = lange_(&NORM, &N, &N, A1, &N, work);
    //BLAS_ge_norm( blas_colmajor, blas_inf_norm, N, N, Residual, N, &Rnorm );
    //BLAS_ge_norm( blas_colmajor, blas_inf_norm, N, N, A1, LDA, &Anorm );

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));

    if ( isnan(Rnorm/(Anorm*N*eps)) || isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}


//--------------------------- check results --------------------

REAL sckres (int n, REAL * A, int lda, REAL * L, int ldl)
{
    REAL zero = 0.0;
    REAL minus_one = -1.0;
    int nminus_one = n - 1;
    REAL one = 1.0;

    REAL nrm, res;
    REAL dummy = 9;

    char NORM = '1';
    char T = 'T';
    char N = 'N';
    char U = 'U';

    res = 0;
    nrm = 0;

    nrm = lange_(&NORM, &n, &n, A, &lda, &dummy);

    laset_(&U, &nminus_one, &nminus_one, &zero, &zero, &L[ldl], &ldl);
    gemm_(&N, &T, &n, &n, &n, &minus_one, L, &ldl, L, &ldl, &one, A, &lda);
    REAL nrm2 = lange_(&NORM, &n, &n, A, &lda, &dummy);

    res = nrm2 / nrm;
    return res;
}


//----------------------------------------------------------------------------------
//			 Changes in data storage
//----------------------------------------------------------------------------------

void print_linear_matrix(int n, REAL *matrix)
{
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%g ", matrix[i*n + j]);
        }
        printf("\n");
    }
}

void read_matrix(char * filename, int n, int ts, REAL *matrix, REAL *checksum)
{
    int i = 0;
    FILE * matrix_file = fopen(filename, "r");
    if (matrix_file != 0) {
    
        while ((i < n*n) && (fscanf(matrix_file, "%g", &matrix[i]) != EOF)) {
            i++;
        }

        // Finished reading matrix; read checksum
        if (fscanf(matrix_file, "%g", checksum) == EOF) {
            printf("Invalid matrix file: could not read checksum\n");
            *checksum = 0.0;
            //exit(1);
        }
    }
    else {
       
        fprintf(stderr, "Matrix file not found, initializing matrix with random values\n");
         
        for (i = 0; i < n*n; i+=n) {
            larnv_(&intONE, &ISEED[0], &n, &matrix[i]);
        }
        
        int j;
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                matrix[j*n + i] = matrix[j*n + i] + matrix[i*n + j];
                matrix[i*n + j] = matrix[j*n + i];
            }
        }
        
        add_to_diag(matrix, n, (REAL) n);
       
        *checksum = 0.0;
    }
     
}
//#include <linux/getcpu.h>
#define _GNU_SOURCE
#include <sched.h>
#include <omp.h>
#if CONVERT_TASK
#ifdef USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task in(N, ts, [N*N]Alin) output([ts*ts]A) label(gather_block) untied
static void gather_block(int N, int ts, REAL *Alin, REAL *A)
{
  //  unsigned cpuid, socketid;
  //  struct getcpu_cache tcache;
    //fprintf(stderr, "CUR CPU: %d\n", sched_getcpu());

   if( CONDITION )
     usleep(DELAY);

   // getcpu(&cpuid, &socketid, &tcache);
   // fprintf(stderr, "CPUID:   %d \n\n", cpuid);
    int i, j;

    for (i = 0; i < ts; i++)
       for (j = 0; j < ts; j++) {
          A[i*ts + j] = Alin[i*N + j];
       }
}
#endif

#if CONVERT_TASK
#ifdef USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task in(ts, N, [ts*ts]A) inout([N*N]Alin) label(scatter_block) untied
static void scatter_block(int N, int ts, REAL *A, REAL *Alin)
{

   if( CONDITION )//|| sched_getcpu() == 2 )//|| sched_getcpu() == 3 || sched_getcpu() == 4)
     usleep(DELAY);

    int i, j;

    for (i = 0; i < ts; i++)
       for (j = 0; j < ts; j++) {
          Alin[i*N + j]= A[i*ts + j];
       }
}
#endif

// static void convert_to_blocks(int ts, int DIM, int N, REAL (*Alin)[N], REAL *(*A)[DIM])
static void convert_to_blocks(int ts, int DIM, int N, REAL Alin[N][N], REAL *A[DIM][DIM])
{
#if CONVERT_TASK
   int i, j;

   for (i = 0; i < DIM; i++)
     for (j = 0; j < DIM; j++) {
         // this for dgemm
         //nanos_current_socket( j % NUM_NODES );
        gather_block ( N, ts, &Alin[i*ts][j*ts], A[i][j]);
     }
#else
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            //A[j/NB][i/NB][(i%NB)*NB+j%NB] = Alin[i*N+j];
            A[j/ts][i/ts][(j%ts)*ts + i%ts] = Alin[j][i];
        }
    }
#endif
}

// static void convert_to_linear(int ts, int DIM, int N, REAL *(*A)[DIM], REAL (*Alin)[N])
static void convert_to_linear(int ts, int DIM, int N, REAL *A[DIM][DIM], REAL Alin[N][N])
{
#if CONVERT_TASK
   int i, j;

   for (i = 0; i < DIM; i++)
     for (j = 0; j < DIM; j++) {
        CHECK_BLOCK_NOT_NULL(A[i][j])
         //nanos_current_socket( j % NUM_NODES );
        scatter_block ( N, ts, A[i][j], (REAL *) &Alin[i*ts][j*ts]);
     }
#else
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Alin[j][i] = A[j/ts][i/ts][(j%ts)*ts + i%ts];
            //Alin[i*N + j] = A[i/NB][j/NB][(j%NB)*NB + i%NB];
        }
    }
#endif
}

#if CONVERT_REQUEST
static REAL * malloc_block (int ts)
{
   REAL *block;
   block = (REAL *) malloc(ts * ts * sizeof(REAL));

   if ( block == NULL ) {
      printf( "ALLOCATION ERROR (Ah block of %d elements )\n", ts );
      exit(-1);
   }

   return block;
}
#endif

//----------------------------------------------------------------------------------
//			 TASKS FOR CHOLESKY
//----------------------------------------------------------------------------------


#if POTRF_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
//this is the case
#pragma omp task inout([NB*NB]A) label(potrf_tile)
#else
#pragma omp task in(NB) inout([NB*NB]A) label(potrf_tile) target device (cuda) untied
#endif
void potrf_tile(REAL *A, int NB)
{

   if( CONDITION )//|| sched_getcpu() == 2 )//|| sched_getcpu() == 3 || sched_getcpu() == 4)
     usleep(DELAY);

    //printf("in portf tile\n");
    char L = 'L';
#if POTRF_SMP
    #if ATLAS
    int INFO;
    int nn=NB;
    spotrf_(&L,
          &nn,
          A,&nn,
          &INFO);
    #else
    mpla
    LONG INFO;
    potrf_(&L, &NB, A, &NB, &INFO);
    #endif
#else
    // Performing Cholesky on GPU
    int INFO;
  #if MAGMA_BLAS
    gpu_potrf(L, NB, A, NB, &INFO);
  #else
    int block = 32;
    REAL * address = A;
    gpu_spotrf_var1_(&L, &NB, (unsigned int *) &address, &NB, &block, &INFO);
  #endif
#endif
}
#if GEMM_SMP
#ifdef USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
//this is the case
#pragma omp task in([NB*NB]A, [NB*NB]B) inout([NB*NB]C) label(gemm_tile) untied
#else
#pragma omp task in([NB*NB]A, [NB*NB]B, NB) inout([NB*NB]C) label(gemm_tile) target device (cuda)
#endif
void gemm_tile(REAL  *A, REAL *B, REAL *C,int NB)
{

   if( CONDITION )// || sched_getcpu() == 2 )//|| sched_getcpu() == 3 || sched_getcpu() == 4)
     usleep(DELAY);

    //printf("in gemm tile\n");
    unsigned char TR = 'T', NT = 'N';
    REAL DONE = 1.0, DMONE = -1.0;
    #if GEMM_SMP
    #if ATLAS
       cblas_sgemm(
        CblasColMajor,
        CblasNoTrans, CblasTrans,
        NB, NB, NB,
        -1.0, A, NB,
              B, NB,
         1.0, C, NB);
    #else
        gemm_(&NT, &TR,
                  &NB, &NB, &NB,
                  &DMONE, A, &NB,
                  B, &NB,
                  &DONE,
                  C, &NB);
     #endif //ATLAS
     #else
     gpu_blas_gemm(NT, TR,
                  NB, NB, NB,
                  DMONE, A, NB,
                  B, NB,
                  DONE,
                  C, NB);
   #endif
}
//---------------------------------------trsm kernal
#if TRSM_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
//this is the case
#pragma omp task in([NB*NB]T) inout([NB*NB]B) label(trsm_tile)
#else
#pragma omp task in([NB*NB]T, NB) inout([NB*NB]B) label(trsm_tile) target device (cuda)
#endif
void trsm_tile(REAL *T, REAL *B, int NB)
{

   if( CONDITION )//|| sched_getcpu() == 2 )//|| sched_getcpu() == 3 || sched_getcpu() == 4)
     usleep(DELAY);

    //printf("in trsm tile!\n");
    char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
    REAL DONE = 1.0;
    #if TRSM_SMP
    #if ATLAS
    cblas_strsm(
       CblasColMajor,
       CblasRight, CblasLower, CblasTrans, CblasNonUnit,
        NB, NB,
        1.0, T, NB,
        B, NB);
    #else //MKL
    trsm_(&RI, &LO, &TR, &NU, &NB, &NB, &DONE, T, &NB, B, &NB );
    #endif //ATLAS
    // Performing STRSM on GPU
    #else
    gpu_blas_trsm(RI, LO, TR, NU, NB, NB,
                  DONE, T, NB, B, NB );
   #endif
}
#if SYRK_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
//this is the case
#pragma omp task in([NB*NB]A) inout([NB*NB]C) label(syrk_tile)
#else
#pragma omp task in([NB*NB]A, NB) inout([NB*NB]C) label(syrk_tile) target device (cuda)
#endif
void syrk_tile(REAL *A, REAL *C, int NB)
{

//   if( CONDITION )//|| sched_getcpu() == 2 )//|| sched_getcpu() == 3 || sched_getcpu() == 4)
//     usleep(DELAY);
//    int x = 0;
//    for(int i = 0; i<5999999999; i++) x = i*x+2;

//    int y = x+x;
    //printf("in syrk tile\n");
    unsigned char LO = 'L', NT = 'N';
    REAL DONE = 1.0, DMONE = -1.0;
    #if POTRF_SMP
    #if ATLAS
    cblas_ssyrk(
        CblasColMajor,
        CblasLower,CblasNoTrans,
        NB, NB,
        -1.0, A, NB,
         1.0, C, NB);
    #else //MKL
    syrk_(&LO, &NT, &NB, &NB, &DMONE, A, &NB, &DONE, C, &NB );
    #endif //ATLAS
    #else
    // Performing SSYRK on GPU
    gpu_blas_syrk(LO, NT, NB, NB,
                DMONE, A, NB, DONE, C, NB );
   #endif
}



//----------------------------------------------------------------------------------
//			 END TASKS FOR CHOLESKY
//----------------------------------------------------------------------------------


#define WARMUP 0 
 
void cholesky(REAL *Alin, REAL** Ah, int ts, int nt)
{
    int i, j, k;
    int* sentinel;
    
#if CONVERT_REQUEST 
    int N = nt * ts;
#endif
#if STOP_SCHED
    NANOS_SAFE( nanos_stop_scheduler() );
    NANOS_SAFE( nanos_wait_until_threads_paused() );
#endif
    //warmup
    #if (WARMUP==1)
      char L = 'L';
      int warmupsize=256;
      #if ATLAS
      int INFO;
      #else
      LONG INFO;
      #endif //ATLAS
    unsigned char  NT = 'N';
    REAL DONE = 1.0, DMONE = -1.0;

      char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
      int warmup2 = 64;

      NT = 'N';
      DMONE = -1.0;
      REAL **W = (REAL**) malloc(sizeof(REAL*)*warmup2);
      for (int i = 0; i < warmup2; i++)  W[i]=(REAL*)malloc(warmup2*sizeof(REAL));

      spotrf_(&L, &warmup2, W[0], &warmup2, &INFO);
 
      L='L';
      //cblas_ssyrk(&L, &NT, &warmup2, &warmup2, &DMONE, W[1], &warmup2, &DONE, W[2], &warmup2 );
      cblas_ssyrk(
        CblasColMajor,
        CblasLower,CblasNoTrans,
        warmup2, warmup2,
        -1.0, W[1], warmup2,
         1.0, W[2], warmup2);


      DONE = 1.0;
      NT='N';
      DMONE = -1.0;
//      cblas_sgemm(&NT, &TR,
//            &warmup2, &warmup2, &warmup2,
//            &DMONE, &W[1], &warmup2,
//            W[0], &warmup2, 
//            &DONE, W[1], &warmup2);
      cblas_sgemm(
        CblasColMajor,
        CblasNoTrans, CblasTrans,
        warmup2, warmup2, warmup2,
        -1.0, W[1], warmup2,
              W[0], warmup2,
         1.0, W[1], warmup2);
//      TR = 'T'; 
//      DONE = -1.0;
//      trsm_( &RI, &LO, &TR, &NU, &warmup2, &warmup2, &DONE, W[0], &warmup2, W[1], &warmup2 );
      fprintf(stderr, "warmup done!\n");    

     /* for (j = 0; j < warmup2; j++) {
         if(W[j]!=NULL)
           free(W[j]);
      }
      free(W);
     */
    #else  
      fprintf(stderr, "No warmup.\n");
    #endif //warmup
    
    // Shuffle across sockets
    for (k = 0; k < nt; k++) {
        TRY_GATHER_BLOCK(N, ts, &Alin[k*ts*N + k*ts], Ah[k*nt + k])
        for (i = 0; i < k; i++) {
            // I don't know why but it works better this way
            //nanos_current_socket( i % NUM_NODES );
            TRY_GATHER_BLOCK(N, ts, &Alin[i*ts*N + i*ts], Ah[i*nt + i])
//#pragma omp task input([ts*ts]Ah[i*nt + k], ts) inout([ts*ts]Ah[k*nt + k]) priority( (nt-i)+10 ) untied
//            CALL_SYRK_TILE (Ah[i*nt + k], Ah[k*nt + k], ts)
            syrk_tile(Ah[i*nt + k], Ah[k*nt + k], ts);
        }
        // Diagonal Block factorization and panel permutations
        //nanos_current_socket( k % NUM_NODES );
//#pragma omp task input(ts) inout([ts*ts]Ah[k*nt + k]) priority( 100000 ) untied
//       CALL_POTRF_TILE(Ah[k*nt + k], ts)
          //task
          potrf_tile(Ah[k*nt + k], ts);
        
        // update trailing matrix
        for (i = k + 1; i < nt; i++) {
            //nanos_current_socket( i % NUM_NODES );
            for (j = 0; j < k; j++) {
                TRY_GATHER_BLOCK(N, ts, &Alin[k*ts*N + j*ts], Ah[k*nt + j])
                TRY_GATHER_BLOCK(N, ts, &Alin[i*ts*N + j*ts], Ah[j*nt + i])
                //CALL_GEMM_TILE(Ah[j*nt + i], Ah[j*nt + k], Ah[k*nt + i], ts)
                //task
                gemm_tile(Ah[j*nt + i], Ah[j*nt + k], Ah[k*nt + i], ts);
            }
            TRY_GATHER_BLOCK(N, ts, &Alin[k*ts*N + i*ts], Ah[k*nt + i])
//#pragma omp task input([ts*ts]Ah[k*nt + k], ts) inout([ts*ts]Ah[k*nt + i]) priority( (nt-i)+10 ) untied
//            CALL_TRSM_TILE (Ah[k*nt + k], Ah[k*nt + i], ts)
              trsm_tile(Ah[k*nt + k], Ah[k*nt + i], ts);
        }
    }
#if STOP_SCHED
    NANOS_SAFE( nanos_start_scheduler() );
    NANOS_SAFE( nanos_wait_until_threads_unpaused() );
#endif
#pragma omp taskwait
}

//--------------------------- MAIN --------------------
int main(int argc, char* argv[])
{
    
    float t1, t2;
    REAL* matrix;
    REAL* original_matrix = NULL;
    REAL *** Ah; 		// Ahdow matrix
    REAL checksum;
    REAL res = -1;
    REAL sum = 0.0;
    int i, j, n, check_result, it;
    int ts, nt, info_factorization;
    char uplo = 'L';
    REAL eps = BLAS_dfpinfo( blas_eps );
    
    struct timeval start;
    struct timeval stop;
    unsigned long elapsed;
    
    if ( argc != 6 ) {
        printf( "cholesky size block_size check_result matrix_file num_iterations\n" );
        exit( -1 );
    }
    n = atoi(argv[1]);
    ts = atoi(argv[2]);
    it = atoi(argv[5]);
    check_result = atoi(argv[3]);
	
    //nanos_get_num_sockets( &NUM_NODES );
    
	// Allocate matrix
    matrix = (REAL *) malloc(n * n * sizeof(REAL));
    if (matrix == NULL) {
        printf("ALLOCATION ERROR\n");
        exit(-1);
    }
    
    read_matrix(argv[4], n, ts, matrix, &checksum);
    
    // Allocate matrix
    if (check_result) {
        original_matrix = (REAL *) malloc(n * n * sizeof(REAL));
        if (original_matrix == NULL) {
            printf("ALLOCATION ERROR\n");
            exit(-1);
        }
    }
    
    nt = n / ts;
    Ah = (REAL***) malloc(it * sizeof(REAL**));
    if(Ah == NULL) {
           printf("ALLOCATION ERROR (Ah)\n");
           exit(-1);
    }
    for(i = 0; i<it; i++)
    {
       // Allocate blocked matrix
       Ah[i] = (REAL **) malloc(nt * nt * sizeof(REAL *));
       if (Ah[i] == NULL) {
           printf("ALLOCATION ERROR (Ah)\n");
           exit(-1);
       }
    }
    for(i = 0; i<it; i++)
    { 
       for (j = 0; j < nt * nt; j++) {
           Ah[i][j]=(REAL *) malloc(ts * ts * sizeof(REAL));
           if (Ah[i][ j ] == NULL) {
               printf("ALLOCATION ERROR (Ah[%d][%d] )\n", i, j);
               exit(-1);
           }
       }
    }
    
    if (check_result) {
        for (i = 0; i < n * n; i++ ) {
            original_matrix[i] = matrix[i];
        }
    }
    
    //int EventSet = PAPI_NULL;
    //unsigned int native = 0x0;
    //if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
//	printf("ERROR in library init!\n");
//	exit(0);
  //  }

    //if(PAPI_create_eventset(&EventSet)!=PAPI_OK) {
//	printf("ERROR creating event set!\n");
//	exit(0);
  //  }
    //if(PAPI_event_name_to_code("rapl:::PP0_ENERGY:PACKAGE0", &native)!=PAPI_OK) {
//	printf("ERROR creating native event\n");
//	exit(0);
  ///  }
    //if(PAPI_add_event(EventSet, native)!=PAPI_OK) {
//	printf("ERROR adding event for Energy!\n");
//	exit(0);	
  //  }
    
    //long_long values[10] = {0};
    //printf("converting to blocks....\n");

#if !CONVERT_REQUEST
    for(i = 0; i<it; i++)
       convert_to_blocks(ts, nt, n, (REAL(*)[n]) matrix, (REAL* (*)[nt]) Ah[i]);
#endif
//    #pragma omp taskwait
    
     t1 = get_time();
//    if(PAPI_start(EventSet) != PAPI_OK) {
//	printf("ERROR in papi start!\n");
//	exit(0);
  //  }


    for(i = 0; i<it; i++)
       cholesky(matrix, Ah[i], ts, nt);

//    #pragma omp taskwait
  //  if (PAPI_stop(EventSet, values) != PAPI_OK) {
//	printf("ERROR in papi stop!\n");
//	exit(0);
  //  }

    t2 = get_time() - t1;
//    if(PAPI_read(EventSet, values)!=PAPI_OK) {
//	printf("ERROR read failed for Energy!\n");
//	exit(0);
//    }


    //printf("Succesfull energy measurement! value = %lld\n", values[0]);
    for(i = 0; i<it; i++)
       convert_to_linear(ts, nt, n, (REAL* (*)[nt]) Ah[i], (REAL (*)[n]) matrix);
    
//    #pragma omp taskwait

 
    if (check_result) {
        info_factorization = check_factorization( n, original_matrix, matrix, n, uplo, eps);
        for (i = 0; i < n * n; i++) {
            sum += original_matrix[i];
        }
        free(original_matrix);
    }
    
    
    float time = t2;
    float gflops = (((1.0 / 3.0) * n * n * n) / ((time) * 1.0e+9));

    // Print configuration
    fprintf(stderr, "\tCONVERT_TASK %d\n\tCONVERT_REQUEST %d\n\tPOTRF_SMP %d\n\tPOTRF_NESTED %d\n\tMAGMA_BLAS %d\n\tTRSM_SMP %d\n\tTRSM_NESTED %d\t\n\tGEMM_SMP %d\n\tGEMM_NESTED %d\n\tSYRK_SMP %d\n\tSYRK_NESTED %d\n\t", CONVERT_TASK, CONVERT_REQUEST, POTRF_SMP, POTRF_NESTED, MAGMA_BLAS,TRSM_SMP,TRSM_NESTED,SYRK_SMP,SYRK_NESTED,GEMM_SMP,GEMM_NESTED);
    
    // Print results
    printf( "============ CHOLESKY RESULTS ============\n" );
    printf( "  matrix size:          %dx%d\n", n, n);
    printf( "  block size:           %dx%d\n", ts, ts);
    printf( "  time (s):             %f\n", time);
    printf( "  performance (gflops): %f\n", gflops);
    printf( "==========================================\n" );

//    printf( "%f ", time);

    
    // Free blocked matrix
    for(i = 0; i<it; i++)
    {
       
       for (j = 0; j < nt * nt; j++) {
          free(Ah[i][j]);
       }
       free(Ah[i]);
    }
    free(Ah);
    free(matrix);
       return 0;
}
