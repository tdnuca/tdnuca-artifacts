#include "qr-tile.hh"
#include <algorithm>

#ifndef _NANOS_H_
inline static void nanos_current_socket(int x) {}
inline static void nanos_get_num_sockets(int *x)
{
  *x = 1;
}
#endif

const char * QR::Left = "L";
const char * QR::Right = "R";
const char * QR::Transpose = "T";
const char * QR::NoTrans= "N";
const char * QR::Forward = "F";
const char * QR::Backward = "B";
const char * QR::Colwise = "C";

QR::QR(Integer M, Integer N, Integer BS)
  : A(M, N, BS),
    T(M, std::min(M, N), BS),
    Q(NULL),
    num_sockets(1),
    BS(BS),
    LAPACK_bs(blocksize(BS))
{
  if (LAPACK_bs <= 1 or LAPACK_bs > BS)
    LAPACK_bs = BS;

  nanos_get_num_sockets(&num_sockets);
  nanos_current_socket(0);
}

QR::~QR()
{
  if (Q != NULL)
    delete Q;
}

void QR::clearMatrix()
{
  for (Integer i = 0; i < A.MB; ++i) {
    nanos_current_socket(i%num_sockets);
    for (Integer j = 0; j < A.NB; ++j) {
      A.blockToZero(i, j);
    }
  }

  // for (Integer i = 0; i < T.MB; ++i) {
  //   nanos_current_socket(i%num_sockets);
  //   for (Integer j = 0; j <= std::min(i,T.NB-1); ++j) {
  //     T.blockToZero(i, j);
  //   }
  // }

  if (Q != NULL) {
    delete Q;
    Q = NULL;
  }
}

void QR::run(bool computeQ)
{
  Integer maxK = std::min(A.NB, A.MB);
  Matrix V(A.BS, A.BS*maxK, A.BS);
  
  for (Integer k = 0; k < maxK; ++k) {
    nanos_current_socket(k%num_sockets);
    this->dgeqrt(A(k, k), V(0, k), T(k, k));
    for (Integer j = k + 1; j < A.NB; ++j) {
      this->dlarfb(V(0, k), T(k, k), A(k, j));
    }
    
    for (Integer i = k + 1; i < A.MB; ++i) {
      nanos_current_socket(i%num_sockets);
      this->dtpqrt(A(k, k), A(i, k), T(i, k));
      for (Integer j = k + 1; j < A.NB; ++j) {
	this->dtpmqrt(A(i, k), T(i, k), A(k, j), A(i, j));
      }
    }
  }

  if (computeQ)
    this->computeQ();

#pragma omp taskwait
}

void QR::computeQ()
{
  if (Q == NULL) {
    Q = new Matrix(A.M, A.M, BS);
  }

  for (Integer i = 0; i < Q->MB; ++i) {
    for (Integer j = 0; j < Q->NB; ++j) {
      nanos_current_socket(j%num_sockets);
      if (i != j)
      	Q->blockToZero(i, j);
      else
	Q->blockToIdentity(i, j);
    }
  }

  Integer maxK = std::min(A.NB, A.MB);
  for (Integer k = 0; k < maxK; ++k) {
    nanos_current_socket(k%num_sockets);
    for (Integer j = 0; j <= k /*Q->NB*/; ++j) {
      this->dlarfb(A(k, k), T(k, k), (*Q)(j, k), true);
    }
    
    for (Integer i = k + 1; i < Q->NB; ++i) {
      nanos_current_socket(i%num_sockets);
      for (Integer j = 0; j < Q->MB; ++j) {
    	this->dtpmqrt(A(i, k), T(i, k), (*Q)(j, k), (*Q)(j, i), true);
      }
    }
  }

  
#pragma omp taskwait
}

void QR::dgeqrt(double *Akk, double *Vkk, double *Tkk)
{
  Integer BS = this->BS;
  Integer LAPACK_bs = this->LAPACK_bs;
#pragma omp task label(dgeqrt) inout([BS*BS]Akk) out([BS*BS]Vkk, [BS*BS]Tkk) firstprivate(BS, LAPACK_bs)
  {
    double *work;
    allocate_aligned(&work, BS*LAPACK_bs);
    Integer info;
    dgeqrt_(&BS, &BS, &LAPACK_bs,
	    Akk, &BS,
	    Tkk, &BS,
	    work,
	    &info);
    free(work);

    if (Akk != Vkk) {
      std::copy(&Akk[0], &Akk[BS*BS], Vkk);
    }
  }
}

void QR::dlarfb(double *Vkk, double *Tkk, double *Akj, bool Q)
{
  Integer BS = this->BS;
  Integer LAPACK_bs = this->LAPACK_bs;
#pragma omp task label(dlarfb) inout([BS*BS]Akj) in([BS*BS]Vkk, [BS*BS]Tkk) firstprivate(BS, LAPACK_bs)
  {
    double *work;
    allocate_aligned(&work, BS*LAPACK_bs);
    dlarfb_(Q ? Right : Left, Q ? NoTrans : Transpose, Q ? Forward : Forward, Colwise,
	    &BS, &BS, &LAPACK_bs,
	    Vkk, &BS,
	    Tkk, &BS,
	    Akj, &BS,
	    work, &BS);
    free(work);
  }
}

void QR::dtpqrt(double *Rkk, double *Aik, double *Tik)
{
  Integer BS = this->BS;
  Integer LAPACK_bs = this->LAPACK_bs;
#pragma omp task label(dtpqrt) inout([BS*BS]Rkk, [BS*BS]Aik, [BS*BS]Tik) firstprivate(BS, LAPACK_bs)
  {
    double *work;
    allocate_aligned(&work, BS*LAPACK_bs);
    Integer l = 0;
    Integer info = 0;
    dtpqrt_(&BS, &BS, &l, &LAPACK_bs,
	    Rkk, &BS,
	    Aik, &BS,
	    Tik, &BS,
	    work,
	    &info);
    free(work);
  }
}

void QR::dtpmqrt(double *Vik, double *Tik, double *Akj, double *Aij, bool Q)
{
  Integer BS = this->BS;
  Integer LAPACK_bs = this->LAPACK_bs;
#pragma omp task label(dtpmqrt) inout([BS*BS]Akj, [BS*BS]Aij) in([BS*BS]Vik, [BS*BS]Tik) firstprivate(BS, LAPACK_bs)
  {
    double *work;
    allocate_aligned(&work, BS*LAPACK_bs);
    Integer l = 0;
    Integer info = 0;
    dtpmqrt_(Q ? Right : Left, Q ? NoTrans : Transpose,
	     &BS, &BS, &BS, &l, &LAPACK_bs,
	     Vik, &BS,
	     Tik, &BS,
	     Akj, &BS,
	     Aij, &BS,
	     work,
	     &info);
    free(work);
  }
}
