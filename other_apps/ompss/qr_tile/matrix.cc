#include "matrix.hh"
#include <iostream>
#include <algorithm>

Matrix::Matrix(long M, long N, long BS)
  : M(M),
    N(N),
    BS(BS),
    MB((M+BS - 1)/BS),
    NB((N+BS - 1)/BS)
{
  allocate(&A, MB*NB);
  if (A == NULL) {
    std::cerr << "Could not allocate matrix" << std::endl;
    exit(1);
  }
  for (long i = 0; i < MB; ++i) {
    for (long j = 0; j < NB; ++j) {
      allocate_aligned(&A[i*NB+j], BS*BS);
      if (A[i*NB+j] == NULL) {
	std::cerr << "Could not allocate tile " << i << ' ' << j << std::endl;
	exit(1);
      }
    }
  }
}

Matrix::~Matrix()
{
  for (long i = 0; i < MB; ++i) {
    for (long j = 0; j < NB; ++j) {
      free(A[i*NB+j]);
    }
  }

  free(A);
}

double *Matrix::operator()(long I, long J)
{
  return A[I*NB+J];
}

void Matrix::blockToZero(long I, long J)
{
#pragma omp task label(zeroBlock) out([BS*BS](A[I*NB+J]))
  std::fill(&A[I*NB+J][0], &A[I*NB+J][BS*BS], 0.0);
}

void Matrix::blockToIdentity(long I, long J)
{
#pragma omp task label(identBlock) out([BS*BS](A[I*NB+J]))
  {
    std::fill(&A[I*NB+J][0], &A[I*NB+J][BS*BS], 0.0);
    for (long i = 0; i < BS; ++i) {
      //std::fill(&A[I*NB+J][i*BS], &A[I*NB+J][i*BS+i], 0.0);
      A[I*NB+J][i*BS+i] = 1.0;
      //std::fill(&A[I*NB+J][i*BS+i+1], &A[I*NB+J][i*BS+BS], 0.0);
    }
  }
}

double &Matrix::elem(long i, long j, bool colwise)
{
  long I = i/BS;
  long J = j/BS;
  long ii = i%BS;
  long jj = j%BS;
  if (colwise) {
    jj *= BS;
  }
  else {
    ii *= BS;
  }

  return A[I*NB+J][ii + jj];
}
