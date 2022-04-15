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
#define CONVERT_TASK 1
#define CONVERT_REQUEST 0
#define POTRF_SMP 1
#define POTRF_NESTED 0
#define TRSM_SMP 1
#define TRSM_NESTED 0
#define SYRK_SMP 1
#define SYRK_NESTED 0
#define GEMM_SMP 1
#define GEMM_NESTED 0
#define TRTRI_SMP 1
#define TRTRI_NESTED 0
#define TRMM_SMP 1
#define TRMM_NESTED 0
#define LAUUM_SMP 1
#define LAUUM_NESTED 0
#define CHOLESKY_INVERSE_SMP 1
#define MAGMA_BLAS 0
#define USE_AFFINITY 1
#define MKL 1

int NUM_NODES;
//#define nanos_current_socket(x) { }

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
#ifdef USE_MKL
#include <mkl.h>
//#else
//#include <cublas.h>
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
#define CALL_TRSM_TILE(x1, x2, x3, x4, x5, x6, x7, priority) omp_trsm(x1, x2, x3, x4, x5, x6, x7);
#else
#define CALL_TRSM_TILE(x1, x2, x3, x4, x5, x6, x7, priority) trsm_tile(x1, x2, x3, x4, x5, x5, x7, priority);
#endif

// Define macro's to make the code cleaner
#if LAUUM_NESTED
#define CALL_LAUUM_TILE(x1, x2) omp_lauum(x1, x2);
#else
#define CALL_LAUUM_TILE(x1, x2) lauum_tile(x1, x2);
#endif

#if SYRK_NESTED
#define CALL_SYRK_TILE(x1, x2, x3, x4, x5, x6, priority) omp_syrk(x1, x2, x3, x4, x5, x6);
#else
#define CALL_SYRK_TILE(x1, x2, x3, x4, x5, x6, priority) syrk_tile(x1, x2, x3, x4, x5, x6, priority);
#endif

#if GEMM_NESTED
#define CALL_GEMM_TILE(x1, x2, x3, x4, x5, x6, x7, x8) omp_gemm(x1, x2, x3, x4, x5, x6, x7, x8);
#else
#define CALL_GEMM_TILE(x1, x2, x3, x4, x5, x6, x7, x8) gemm_tile(x1, x2, x3, x4, x5, x6, x7, x8);
#endif

#if TRTRI_NESTED
#define CALL_TRTRI_TILE(x1, x2, x3) omp_dtrtri(x1, x2, x3);
#else
#define CALL_TRTRI_TILE(x1, x2, x3) trtri_tile(x1, x2, x3);
#endif

#if TRMM_NESTED
#define CALL_TRMM_TILE(x1, x2, x3, x4, x5, x6, x7, priority) omp_trmm(x1, x2, x3, x4, x5, x6, x7);
#else
#define CALL_TRMM_TILE(x1, x2, x3, x4, x5, x6, x7, priority) trmm_tile(x1, x2, x3, x4, x5, x6, x7, priority);
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
#define trtri_      dtrtri_
#define lauum_      dlauum_
#define symm_       dsymm_
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
#define  trtri_      strtri_
#define  lauum_      slauum_
#define  symm_       ssymm_
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
void laset_(const char * UPLO, const int * M, const int * N, const REAL * ALPHA, const REAL * BETA, REAL * A, const int * LDA);

void gemm_ (const char *transa, const char *transb, int *l, int *n, int *m, REAL *alpha,
    const void *a, int *lda, void *b, int *ldb, REAL *beta, void *c, int *ldc);

REAL lange_ (const char *norm, const int *m, const int *n, const REAL *a, const int *lda, REAL *work);
void lacpy_ (const char *norm, const int *m, const int *n, const REAL *a, const int *lda, REAL *b, const int *ldb);
void trtri_(const char *uplo, const char *diag, const int *n, REAL *a, const int *lda, LONG *info);

void larnv_ (const int *idist, int *iseed, const int *n, REAL *x);


void potrf_(const char* uplo, const int* n, REAL* a, const int* lda, LONG* info );
void trsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, REAL *alpha, REAL *a, int *lda, REAL *b, int *ldb);
void trmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, REAL *alpha, REAL *a, int *lda, REAL *b, int *ldb);
void syrk_(char *uplo, char *trans, int *n, int *k, REAL *alpha, REAL *a, int *lda, REAL *beta, REAL *c, int *ldc);
void symm_(const char *side, const char *uplo, const int *m, const int *n, REAL *alpha, REAL *a, int *lda, REAL*b, int *ldb, REAL *beta, REAL *c, int *ldc);

void lauum_(const char *uplo, const int *n, REAL *a, const int *lda, int *info);

float get_time();
static int check_factorization(int, REAL*, REAL*, int, char , REAL);
static int check_inverse(int, REAL *, REAL *, int, char, REAL);

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

void cholesky_inverse(REAL *Alin, REAL** Ah, int ts, int nt, int matrix_size);

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




//-----------------------check inverse---------------------------

static int check_inverse(int N, REAL *A1, REAL *A2, int LDA, char uplo, REAL eps )
{
  int info_inverse;
  int i, j;
  REAL Rnorm, Anorm, Ainvnorm, result;
  REAL alpha, beta, zone;
  REAL *work = (REAL *)malloc(N*N*sizeof(REAL));

  char NORM = 'I', ALL = 'A', UP = 'U', LO = 'L', TR = 'T', NU = 'N', RI = 'R';
  alpha = -1.0;
  beta  = 0.0;
  zone = 1.0;

  /* Rebuild the other part of the inverse matrix */
  if (uplo == 'U'){
    for(j=0; j<N; j++)
      for(i=0; i<j; i++)
        *(A2+j+i*LDA) = *(A2+i+j*LDA);
    symm_(&LO, &uplo, &N, &N, &alpha, A2, &LDA, A1, &LDA, &beta, work, &N);

  }
  else {
    for(j=0; j<N; j++)
      for(i=j; i<N; i++)
        *(A2+j+i*LDA) = *(A2+i+j*LDA);
    symm_(&RI, &uplo, &N, &N, &alpha, A2, &LDA, A1, &LDA, &beta, work, &N);
  }

  /* Add the identity matrix to work */
  for(i=0; i<N; i++)
    *(work+i+i*N) = *(work+i+i*N) + zone;

  Rnorm = lange_(&NORM, &N, &N, work, &N, work);
  Anorm = lange_(&NORM, &N, &N, A1, &N, work);
  Ainvnorm = lange_(&NORM, &N, &N, A2, &N, work);

  // BLAS_dge_norm( blas_colmajor, blas_one_norm, N, N, work, N, &Rnorm );
  //BLAS_dge_norm( blas_colmajor, blas_one_norm, N, N, A1, LDA, &Anorm );
  //BLAS_dge_norm( blas_colmajor, blas_one_norm, N, N, A2, LDA, &Ainvnorm );

  if (getenv("PLASMA_TESTING_VERBOSE"))
    printf( "||A||_1=%f\n||Ainv||_1=%f\n||Id - A*Ainv||_1=%e\n", Anorm, Ainvnorm, Rnorm );

  result = Rnorm / ( (Anorm*Ainvnorm)*N*eps ) ;
  printf("============\n");
  printf("Checking the Residual of the inverse \n");
  printf("-- ||Id - A*Ainv||_1/((||A||_1||Ainv||_1).N.eps) = %e \n", result);

  if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
    printf("-- The inverse is suspicious ! \n");
    info_inverse = 1;
  }
  else{
    printf("-- The inverse is CORRECT ! \n");
    info_inverse = 0;
  }

  free(work);

  return info_inverse;
}

// Numa facilities
void* cholesky_malloc( size_t size, unsigned node );
void cholesky_free( void* address, size_t size );
int get_node( void* address );

void* cholesky_malloc( size_t size, unsigned node )
{
  //#ifdef USE_NUMA
  //    return numa_alloc_onnode( size, node );
  //#else
  void* address;
  posix_memalign( &address, getpagesize(), size );
  return address;
  //return malloc( size );
  //#endif
}

void cholesky_free( void* address, size_t size )
{
#ifdef USE_NUMA
  numa_free( address, size );
#else
  free( address );
#endif
}

#ifdef USE_NUMA
int get_node( void* address )
{
  //printf( "mask: %p\n", /*~( 4096ul - 0x1 )*/-4096l);
  //void* page = ( void* ) ((uintptr_t)address & ~( 4096ul - 0x1 ) );
  void* page = ( void* ) ((uintptr_t)address & -4096l );
  //printf( "Address %p was converted to %p\n", address, page );
  /*here you should align ptr_to_check to page boundary */
  int status[1];
  int ret_code;
  status[0]=-1;
  ret_code=numa_move_pages(0 /*self memory */, 1, &page,
      NULL/*nodes*/, status, 0 /*flags*/);
  return status[0];
}
#endif

void check_affinity( int node, void* param1, void* param2, void* param3, const char* func )
{
  return;
#ifdef USE_NUMA
  int node1 = get_node( param1 );
  int node2 = get_node( param2 );
  int node3 = (param3!=NULL) ? get_node( param3 ) : -10 ;

  //fprintf( stderr, "%s: memory at nodes %d, %d, %d, I will use node %d\n", func, node1, node2, node3, node );
  int candidate = node;
  if ( ( node1 == node2 ) || ( node1 == node3 ) )
    candidate = node1;
  else if ( node2 == node3 )
    candidate = node2;
  if ( candidate != node )
    fprintf( stderr, "%s: memory at nodes %d, %d, %d, I will use node %d, should be %d\n", func, node1, node2, node3, node, candidate );
  nanos_current_socket( candidate );
#else
  nanos_current_socket( node );
#endif
}

//-----------------------------------------------

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

    printf("Matrix file not found, initializing matrix with random values\n");

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

#if CONVERT_TASK
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task in([N*(ts-1)+ts]Alin) out([ts*ts]A) 
static void gather_block(int N, int ts, REAL *Alin, REAL *A)
{
  int i, j;

  for (i = 0; i < ts; i++)
    for (j = 0; j < ts; j++) {
      A[i*ts + j] = Alin[i*N + j];
    }
}
#endif

#if CONVERT_TASK
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task in([ts*ts]A) out([N*(ts-1)+ts]Alin)
static void scatter_block(int N, int ts, REAL *A, REAL *Alin)
{
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
      nanos_current_socket( j % NUM_NODES );
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
        nanos_current_socket( j % NUM_NODES );
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
#pragma omp task firstprivate(NB) inout([NB*NB]A) priority(100000) untied
#else
#pragma omp task firstprivate(NB) inout([NB*NB]A) target device (cuda)
#endif
void potrf_tile(REAL *A, int NB)
{
  char L = 'L';
#if POTRF_SMP
  LONG INFO;
  potrf_(&L, &NB, A, &NB, &INFO);
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

//--------------------------------------

#if LAUUM_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(NB) inout([NB*NB]A) untied priority( 100000 )
#else
#pragma omp task firstprivate(NB) inout([NB*NB]A) target device (cuda)
#endif
void lauum_tile(REAL *A, int NB)
{
  char L = 'L';
  LONG INFO;
  lauum_(&L, &NB, A, &NB, &INFO);

}

//-------------------------------------

#if TRMM_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(NB) in([NB*NB]T) inout([NB*NB]B) priority( priority ) untied
#else
#pragma omp task firstprivate(NB) in([NB*NB]T) inout([NB*NB]B) target device (cuda)
#endif
void trmm_tile(REAL *T, REAL *B, int NB, char side ,char trans, char diag, REAL alpha, unsigned priority)
{
  char uplo='L';
  trmm_(&side, &uplo, &trans, &diag, &NB, &NB, &alpha, T, &NB, B, &NB);

}

//------------------------------------

#if TRTRI_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(NB) inout([NB*NB]A) priority( 100000 ) untied
#else
#pragma omp task firstprivate(NB) inout([NB*NB]A) target device (cuda)
#endif
void trtri_tile(REAL *A, int NB, char diag)
{
  char L = 'L';
  LONG INFO;
  trtri_(&L, &diag, &NB, A, &NB, &INFO);

}
//-----------------------------------------
#if GEMM_SMP
#ifdef USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(NB) in([NB*NB]A, [NB*NB]B) inout([NB*NB]C) untied
#else
#pragma omp task firstprivate(NB) in([NB*NB]A, [NB*NB]B) inout([NB*NB]C) target device (cuda)
#endif
void gemm_tile(REAL  *A, REAL *B, REAL *C,int NB, char transa, char transb, REAL alpha, REAL beta)
{
#if 1
  nanos_event_t e;
  e.value = NB*NB*NB*2;
  nanos_instrument_register_key ( &e.key, "flops", "Flops", false );
  e.type = NANOS_BURST_START;
  nanos_instrument_events( 1, &e);
#endif
  REAL DONE = 1.0, DMONE = -1.0;
#if GEMM_SMP
  gemm_(&transa, &transb,
      &NB, &NB, &NB,
      &alpha, A, &NB,
      B, &NB,
      &beta,
      C, &NB);
#else
  gpu_blas_gemm(transa, transb,
      NB, NB, NB,
      DMONE, A, NB,
      B, NB,
      DONE,
      C, NB);
#endif
#if 1
  e.type = NANOS_BURST_END;
  nanos_instrument_events( 1, &e);
#endif
}
//---------------------------------------trsm kernal
#if TRSM_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(NB) in([NB*NB]T) inout([NB*NB]B) priority(priority) untied
#else
#pragma omp task firstprivate(NB) in([NB*NB]T) inout([NB*NB]B) priority(priority) target device (cuda)
#endif
void trsm_tile(REAL *T, REAL *B, int NB, char side , char trans, char diag, REAL alpha, unsigned priority)
{
  char LO = 'L';
  REAL DONE = 1.0;
  char di='N';
#if TRSM_SMP
  trsm_(&side, &LO, &trans, &di , &NB, &NB, &alpha, T, &NB, B, &NB);
  // Performing STRSM on GPU
#else
  gpu_blas_trsm(RI, LO, trans, NU, NB, NB,
      DONE, T, NB, B, NB );
#endif
}

//------------------------------------------------

#if SYRK_SMP
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(NB) in([NB*NB]A) inout([NB*NB]C) priority(priority) untied
#else
#pragma omp task firstprivate(NB) in([NB*NB]A) inout([NB*NB]C) priority(priority) target device (cuda)
#endif
void syrk_tile(REAL *A, REAL *C, int NB, char trans, REAL alpha, REAL beta, unsigned priority)
{
  unsigned char LO = 'L';
  REAL DONE = 1.0, DMONE = -1.0;
#if POTRF_SMP
  syrk_(&LO, &trans, &NB, &NB, &alpha, A, &NB, &beta, C, &NB );
#else
  // Performing SSYRK on GPU
  gpu_blas_syrk(LO, NT, NB, NB,
      DMONE, A, NB, DONE, C, NB );
#endif
}



//----------------------------------------------------------------------------------
//			 END TASKS FOR CHOLESKY
//----------------------------------------------------------------------------------

//*****************
//Stage 1: Cholesky factorization
#if 0
#if CHOLESKY_INVERSE_SMP 
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task firstprivate(nt) inout([nt*nt]Ah)
#else
#pragma omp task firstprivate(nt) inout([nt*nt]Ah) target device (cuda)
#endif
#endif
void cholesky_inverse(REAL *Alin, REAL** Ah, int ts, int nt, int matrix_size)
{
  int m, n, k;
  int* sentinel;

#if CONVERT_REQUEST
  int N = nt * ts;
#endif
  // Shuffle across sockets
  for (k = 0; k < nt; k++) {
    nanos_current_socket( k % NUM_NODES );
    CALL_POTRF_TILE(Ah[k*nt + k], ts)
      // Triangular systems
      for ( m= k + 1; m < nt; m++) {
        CALL_TRSM_TILE (Ah[k*nt + k], Ah[k*nt + m], ts, 'R' ,'T', 'N', 1.0, (nt-m)+10)
      }
    // update trailing matrix
    for (m= k + 1; m < nt; m++) {
      nanos_current_socket( m %  NUM_NODES );
      // gmiranda: in the normal cholesky we used to do this
      for (n = k + 1; n < m; n++) {
        check_affinity( m %  NUM_NODES, Ah[k*nt + m], Ah[k*nt + n], Ah[n*nt + m], "gemm_stage1" );
        CALL_GEMM_TILE(Ah[k*nt + m], Ah[k*nt + n], Ah[n*nt + m], ts, 'N', 'T', -1.0, 1.0)
      }
      CALL_SYRK_TILE (Ah[k*nt + m], Ah[m*nt + m], ts, 'N', -1.0, 1.0, (nt-m)+10)
        //            CALL_SYRK_TILE (Ah[k*nt + m], Ah[m*nt + m], ts, 'N', -1.0, 1.0,(nt-m)+10)
        //           for (n = k + 1; n < m; n++) {
        //                CALL_GEMM_TILE(Ah[k*nt + m], Ah[k*nt + n], Ah[n*nt + m], ts, 'N', 'T', -1.0, 1.0)
        //            }

    } // End of stage 1

    // Start of stage 2
    for (m= k + 1; m< nt; m++) {
      nanos_current_socket(m% NUM_NODES);
      CALL_TRSM_TILE (Ah[k*nt + k], Ah[k*nt + m], ts, 'R' ,'N', 'N', -1.0, (nt-m)+10 )
        // gmiranda: doesn't this need to go in another loop?
        for (n = 0; n < k; n++) {
          check_affinity( m %  NUM_NODES, Ah[k*nt + m], Ah[n*nt + k], Ah[n*nt + m], "gemm_stage2" );
          CALL_GEMM_TILE(Ah[k*nt + m], Ah[n*nt + k], Ah[n*nt + m], ts, 'N', 'N', 1.0, 1.0)
        }
    }

    for (m =0; m < k; m++) {
      nanos_current_socket(m % NUM_NODES);
      CALL_TRSM_TILE (Ah[k*nt + k], Ah[m*nt + k], ts, 'L' ,'N', 'N', 1.0, (nt-m)+10)
    }

    nanos_current_socket( k % NUM_NODES );
    CALL_TRTRI_TILE(Ah[k*nt + k], ts, 'N')
      // End of stage 2


      // Start of stage 3
      for (n = 0; n < k; n++) { 
        nanos_current_socket( n % NUM_NODES);
        CALL_SYRK_TILE (Ah[n*nt + k], Ah[n*nt + n], ts, 'T', 1.0, 1.0, (nt-n)+10)
          for (m= n+1; m < k; m++) {
            nanos_current_socket(k % NUM_NODES);
            check_affinity( k %  NUM_NODES, Ah[m*nt + k], Ah[n*nt + k], Ah[n*nt + m], "gemm_stage3" );
            CALL_GEMM_TILE(Ah[m*nt + k], Ah[n*nt + k], Ah[n*nt + m], ts, 'T', 'N', 1.0, 1.0)
          }
      }

    for (n =0; n < k; n++) {
      nanos_current_socket(n % NUM_NODES);
      CALL_TRMM_TILE (Ah[k*nt + k], Ah[n*nt + k], ts, 'L','T', 'N', 1.0, (nt-n)+10)
    }

    nanos_current_socket(k % NUM_NODES);
    CALL_LAUUM_TILE(Ah[k*nt + k], ts)
      // End of stage 3


  }


  //#pragma omp taskwait noflush
}

// MKL and others need to allocate buffers before computing.
void warmup()
{
  int block = 512;
  REAL *A = (REAL*) malloc( block * block * sizeof( REAL ) );
  REAL *B = (REAL*) malloc( block * block * sizeof( REAL ) );
  REAL *C = (REAL*) malloc( block * block * sizeof( REAL ) );

  char transa = 'N';
  char transb = 'N';

  REAL DONE = 1.0, DMONE = -1.0;

  gemm_(&transa, &transb,
      &block, &block, &block,
      &DONE, A, &block,
      B, &block,
      &DMONE,
      C, &block);
  free( A );
  free( B );
  free( C );

}

//--------------------------- MAIN --------------------
int main(int argc, char* argv[])
{

  float t1, t2;
  REAL* matrix;
  REAL* original_matrix = NULL;
  REAL ** Ah; 		// Ahdow matrix
  REAL checksum;
  REAL res = -1;
  REAL sum = 0.0;
  int i, j, check_result;
  size_t n;
  int ts, nt, info_factorization, info_inverse;
  char uplo = 'L';
  REAL eps = BLAS_dfpinfo( blas_eps );
  bool stealing;


  if ( argc != 5 ) {
    printf( "Usage: inverse size block_size check_result matrix_file\n" );
    exit( -1 );
  }
  n = atoi(argv[1]);
  ts = atoi(argv[2]);
  check_result = atoi(argv[3]);
  if((n%ts)!=0) exit(-1);

  nanos_get_num_sockets( &NUM_NODES );
  fprintf( stderr, "Running with %d nodes\n", NUM_NODES );

  // Allocate matrix
  matrix = (REAL *) malloc(n * n * sizeof(REAL));
  if (matrix == NULL) {
    printf("ALLOCATION ERROR\n");
    exit(-1);
  }
#if NANOS_API_COPIES_API >= 1004
#pragma omp register( [n*n]matrix )
#endif

  warmup();

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

  // Allocate blocked matrix
  Ah = (REAL **) malloc(nt * nt * sizeof(REAL *));
  if (Ah == NULL) {
    printf("ALLOCATION ERROR (Ah)\n");
    exit(-1);
  }

  for (j = 0; j < nt * nt; j++) {
    Ah[j]=(REAL *) cholesky_malloc(ts * ts * sizeof(REAL), j%NUM_NODES);
    if (Ah[ j ] == NULL) {
      printf("ALLOCATION ERROR (Ah[%d] )\n", j);
      exit(-1);
    }
    REAL * block = Ah[j];
#if NANOS_API_COPIES_API >= 1004
#pragma omp register( [ts*ts]block )
#endif
  }

  if (check_result) {
    for (i = 0; i < n * n; i++ ) {
      original_matrix[i] = matrix[i];
    }
  }
  if (check_result)
  {
    nanos_scheduler_get_stealing( &stealing );
    // Disable stealing so that the memory is evenly distributed
    if ( stealing ) nanos_scheduler_set_stealing( false );
#if !CONVERT_REQUEST

    convert_to_blocks(ts, nt, n, (REAL(*)[n]) matrix, (REAL* (*)[nt]) Ah);
#endif
#pragma omp taskwait noflush
    // Enable it again
    if ( stealing ) nanos_scheduler_set_stealing( true );

    cholesky_inverse(matrix, Ah, ts, nt, n);

    //#pragma omp taskwait noflush

    //convert_to_linear(ts, nt, n, (REAL* (*)[nt]) Ah, (REAL (*)[n]) matrix);

    //#pragma omp taskwait noflush

    //info_factorization = check_factorization( n, original_matrix, matrix, n, uplo, eps);

    //-------------------------------

    //------check inverse------------


#pragma omp taskwait noflush

    convert_to_linear(ts, nt, n, (REAL* (*)[nt]) Ah, (REAL (*)[n]) matrix);

#pragma omp taskwait noflush

    info_inverse = check_inverse(n, original_matrix, matrix, n, uplo, eps);
    //-----------------------------------------
    free(original_matrix);

  }
  else
  {
    nanos_scheduler_get_stealing( &stealing );
    // Disable stealing so that the memory is evenly distributed
    if ( stealing ) nanos_scheduler_set_stealing( false );

#if !CONVERT_REQUEST
    convert_to_blocks(ts, nt, n, (REAL(*)[n]) matrix, (REAL* (*)[nt]) Ah);
#endif
#pragma omp taskwait noflush
    // Enable it again
    if ( stealing ) nanos_scheduler_set_stealing( true );

    t1 = get_time();

    cholesky_inverse(matrix, Ah, ts, nt, n);

#pragma omp taskwait noflush

    t2 = get_time() - t1;

    convert_to_linear(ts, nt, n, (REAL* (*)[nt]) Ah, (REAL (*)[n]) matrix);

#pragma omp taskwait noflush


    float time = t2;
    /*
     * note (gmiranda): lawn41 expects spotri to receive a factorised A matrix.
     * As we compute the Cholesky factorization, we must take that into account.
     */
    float gflops = (1.0 / 3.0) * n * n * n + ( 1.0 / 2.0 ) * n * n + ( 1.0 / 6.0 ) * n;
    gflops += (((2.0 / 3.0) * n * n * n)+((1.0/2.0)*n*n)+((5.0/6.0)*n));
    gflops /= ((time) * 1.0e+9);

    // Print configuration
    //printf("\tCONVERT_TASK %d\n\tCONVERT_REQUEST %d\n\tPOTRF_SMP %d\n\tPOTRF_NESTED %d\n\tMAGMA_BLAS %d\n\tTRSM_SMP %d\n\tTRSM_NESTED %d\t\n\tGEMM_SMP %d\n\tGEMM_NESTED %d\n\tSYRK_SMP %d\n\tSYRK_NESTED %d\n\t", CONVERT_TASK, CONVERT_REQUEST, POTRF_SMP, POTRF_NESTED, MAGMA_BLAS,TRSM_SMP,TRSM_NESTED,SYRK_SMP,SYRK_NESTED,GEMM_SMP,GEMM_NESTED);

    // Print results
    printf( "============ CHOLESKY RESULTS ============\n" );
    printf( "  matrix size:          %dx%d\n", n, n);
    printf( "  block size:           %dx%d\n", ts, ts);
    printf( "  time (s):             %f\n", time);
    printf( "  performance (gflops): %f\n", gflops);
    printf( "==========================================\n" );
  }    

  // Free blocked matrix
  for (j = 0; j < nt * nt; j++) {
    cholesky_free(Ah[j], ts * ts * sizeof(REAL));
  }
  free(Ah);
  free(matrix);
  return 0;
}
