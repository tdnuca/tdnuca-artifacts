#ifndef __QR_HH
#define __QR_HH
#include "matrix.hh"
#include "lapack_headers.h"

struct QR
{
  Matrix A;
  Matrix T;
  Matrix *Q;

  int num_sockets;
  Integer BS;
  Integer LAPACK_bs;

  static const char * Left;
  static const char * Right;
  static const char * Transpose;
  static const char * NoTrans;
  static const char * Forward;
  static const char * Backward;
  static const char * Colwise;

  QR(Integer M, Integer N, Integer BS);
  ~QR();
  void clearMatrix();
  void run(bool computeQ=false);
  void computeQ();

private:

  void dgeqrt(double *Akk, double *Vkk, double *Tkk);
  void dlarfb(double *Vkk, double *Tkk, double *Akj, bool Q=false);
  void dtpqrt(double *Rkk, double *Aik, double *Tik);
  void dtpmqrt(double *Vik, double *Tik, double *Akj, double *Aij, bool Q=false);
}; 

#endif
