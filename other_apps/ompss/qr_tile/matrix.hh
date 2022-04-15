#ifndef __MATRIX_HH
#define __MATRIX_HH
#include "allocation.hh"

struct Matrix
{
  long M;
  long N;
  long BS;
  long MB;
  long NB;

  char concHelp;

  double **A; // A matrix
  Matrix(long M, long N, long BS);

  double *operator()(long I, long J);

  void blockToZero(long I, long J);

  void blockToIdentity(long I, long J);

  double &elem(long i, long j, bool colwise=true);

  virtual ~Matrix();
};

#endif
