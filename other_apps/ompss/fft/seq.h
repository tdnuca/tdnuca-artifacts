/********************************************************************
   F A S T   F O U R I E R   T R A N S F O R M   P R O G R A M S

  by Wang Jian-Sheng 4 Nov 1998, added fft2D(), 11 Apr 2003
---------------------------------------------------------------------

  Reference: "Computational Frameworks for the Fast Fourier
              Transform", Charles Van Loan, SIAM, 1992.

  There are many FFT algorithms, the most important ones are
     COOLEY-TUKEY:  in place, bit reversal
     STOCKHAM AUTOSORT:  additional memory size of input data
     MIXED RADIX:  20% less operations comparing to Cooley-Tukey
     PRIME FACTOR: arbitrary length n

  We use a combination of the Stockham autosort algorithm 1.7.2,
  page 57, and multirow Cooley-Tukey (3.1.7), page 124, of the
  reference above.

  The discrete Fourier transform is defined by
  y[k] = sum_(j=0,n-1) x[j] exp(-2 Pi sqrt[-1] j k/n),
  k=0,1,...,n-1.  The factor (1/n) is not included.
  If y[]<-x[]; fft(x,n,1); fft(x,n,-1); then y[]==x[]/n is true.
  Three dimensional transform is generalized straightforwardly.

   Interface and usage:
   1D Fourier transform
   Use: fft(x, n, flag)
      x    : an array of structure type complex;
      n    : the size of data, must be a power of 2;
      flag : 1 for forward transform, -1 for inverse transform.

   3D Fourier transform
   Use :  fft3D(x, n1, n2, n3, flag)
     x    : 1D array of type complex representing 3D array;
            mapping through C convention, i.e.,
            (i,j,k) -> k + n3*j + n2*n3*i;
     n1, n2, n3 : dimensions in three directions;
     flag : same as in 1D.

   2D FFT is similar but with n1 and n2 only.

**********************************************************************/

/* Inclusion of standard C libraries */

#ifndef REF_H_
#define REF_H_

#include <math.h>
#include <complex.h>


#define  REALSIZE     8                                  /* in units of byte */
typedef double real;                 /* can be long double, double, or float */


/* Mathematical functions and constants */
#undef M_PI
#if (REALSIZE==16)
#define sin  sinl
#define cos  cosl
#define fabs  fabsl
#define M_PI 3.1415926535897932384626433832795L
#else
#define M_PI 3.1415926535897932385E0
#endif

void seq_fft(int n, complex double *x, int flag);
void seq_fft2D(int n1, int n2, complex double (*x)[n2], int flag);
void seq_fft3D(int n1, int n2, int n3, complex double (*x)[n2][n3], int flag);

#endif //
