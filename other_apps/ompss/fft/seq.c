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

#include <math.h>
#include <complex.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <err.h>

#include "fft.h"


static inline void swap_ptr(void** x, void** y)
{
	void *t = *x;
	*x = *y;
	*y = t;
}


static inline void swap(complex double* x, complex double* y)
{
	complex double t = *x;
	*x = *y;
	*y = t;
}


static inline complex double * compute_omegas(const int n, int flag)
{
	complex double *w = malloc(n * sizeof(*w));

	// precompute the omegas (j, n). Later use omega(j, Ls / 2) = omega(2 * j, Ls) for omegas with smaller n values
	for (int j = 0; j < n; j++)
		w[j] = cos(M_PI * j / n) - I * flag * sin(M_PI * j / n);

	return w;
}

/*----------------------------------------------------------------------*/
/* Truncated Stockham algorithm for single-column vector
       X(n) <- F_n X(n)
   x input data of size n and w the omegas for n / 2.  This function is
   supposed to be internal (static), not used by application.  Note that the
   terminology of column or row respect to algorithms in the Loan's book is
   reversed, because we use row major convention of C.
*/
static void stockham_1D(const int n1, complex double *x, complex double *y, complex double *w)
{
	complex double *x_orig = x;

	// Find a good chunk size, e.g. n1 / 8 for 8 threads.
	/* loops log2(n1) times. Ls=L/2 is the L star. Ls * r == n1 / 2 */
	for (int Ls = 1, r = n1 / 2; r >= 1; r /= 2, Ls *= 2)
	{
		swap_ptr((void**)&x, (void**)&y);

		// Iterate over 0..n1/2 of dimension 1, in chunks of size r and fill x linearly and at pos shited by n1/2 simultaneously
		for (int j = 0; j < Ls; j++)
		{
			int j1 = j * r;
			int j2 = j * r * 2;

			// get y[j2+r..j2+2r] and y[j2..j2+r] into x[j1..j1+r] and x[j1+n1/2..j1+n1/2+r]
			for (int k = 0; k < r; k++)          /* "butterfly" operation */
			{
				complex double t = w[j * r] * y[j2 + k + r];
				x[j1 + k] = y[j2 + k] + t;
				x[j1 + k + n1 / 2] = y[j2 + k] - t;
			}
		}
	}

	/* copy back to permanent memory if not already there: performed iff log2(n1) is odd */
	if (x != x_orig)
	{
		swap_ptr((void**)&x, (void**)&y);
		memcpy(x, y, n1 * sizeof(*y));
	}
}

/*----------------------------------------------------------------------*/
/* Truncated Stockham algorithm for multi-column vector,
       X(n1,n2) <- F_n1 X(n1,n2)
   x input data is a n1 by n2 dimensional array, w are the omegas for n2 / 2.
   This function is supposed to be internal (static), not used by application.
   Note that the terminology of column or row respect to algorithms in the
   Loan's book is reversed, because we use row major convention of C.
*/
static void stockham(const int n1, const int n2, complex double (*x)[n2], complex double (*y)[n2], complex double *w)
{
	complex double (*x_orig)[n2] = x;

	for (int Ls = 1, r = n1 / 2; r >= 1; r /= 2, Ls *= 2)  /* loops log2(n1) times. Ls=L/2 is the L star  */
	{
		swap_ptr((void**)&x, (void**)&y);

		// Ls * r == n1 / 2
		// Iterate over 0..n1/2 of dimension 1, in chunks of size r and fill x linearly and at pos shited by n1/2 simultaneously
		for (int j = 0; j < Ls; j++)
		{
			for (int k = 0; k < r; k++)          /* "butterfly" operation */
			{
				int k1 = j * r     + k;
				int k2 = j * r * 2 + k;
				for (int l = 0; l < n2; l++)
				{
					// get y[j2r+r][0..n2] and y[j2r+r+k][0..n2] into x[jr+k][0..n2] and x[n1/2+jr+k][0..n2]
					complex double t = w[j * r] * y[k2 + r][l];
					x[k1][l] = y[k2][l] + t;
					x[k1 + n1 / 2][l] = y[k2][l] - t;
				}
			}
		}
	}

	/* copy back to permanent memory if not already there: performed iff log2(n1) is odd */
	if (x != x_orig)
	{
		swap_ptr((void**)&x, (void**)&y);
		memcpy(x, y, n1 * n2 * sizeof(x[0][0]));
	}
}


static inline int bit_reversal(int n1, int k)
{
	int j = 0, m = k, p = 1;              /* p = 2^q,  q used in the book */
	while (p < n1)
	{
		j = 2 * j + (m & 1);
		m >>= 1;
		p <<= 1;
	}
	assert(p == n1);                   /* make sure n1 is a power of two */
	return j;
}


/* The Cooley-Tukey multiple column algorithm, see page 124 of Loan.
   x[] is input data, overwritten by output, and is a n1 by n2 array.
   w contains the omegas for n1 / 2.
*/
static inline void cooley_tukey(const int n1, const int n2, complex double (*x)[n2], complex double *w)
{
	/* do bit reversal permutation */
	for (int k = 0; k < n1; k++)         /* This is algorithms 1.5.1 and 1.5.2. */
	{
		int j = bit_reversal(n1, k);
		if (j > k)
			for (int i = 0; i < n2; i++)                      /* swap k <-> j row for all columns */
				swap(x[k] + i, x[j] + i);
	}

	/* This is (3.1.7), page 124 */
	for (int Ls = 1, r = n1 / 2; Ls < n1; Ls *= 2, r /= 2) // Ls * r == n1 / 2
	{
		for (int j = 0; j < Ls; j++)
		{
			// update x[j..j+n1][0..n2] and x[(j..j+n1) + Ls][0..n2]
			for (int k = 0; k < r; k++)
			{
				int j1 = j + Ls * (2 * k);
				int j2 = j + Ls * (2 * k + 1);
				// update x[k][0..n2] and x[k + Ls][0..n2]
				for (int i = 0; i < n2; i++)                       /* for each row */
				{
					complex double t = w[j * r] * x[j2][i];
					x[j2][i] = x[j1][i] - t;
					x[j1][i] += t;
				}
			}
		}
	}
}



/* 1D Fourier transform:
   Simply call stockham with proper arguments.
   Allocated working space of size n dynamically.
*/
void seq_fft(int n, complex double *x, int flag)
{
	complex double *y = malloc(n * sizeof(*y));
	if (!y) err(1, "Failed to allocate y for %d complex doubles", n);

	complex double *w = compute_omegas(n / 2, flag);

	stockham_1D(n, x, y, w);

	free(y);
	free(w);
}


void seq_fft2D(int n1, int n2, complex double (*x)[n2], int flag)
{
	complex double (*y)[n2] = (complex double (*)[n2])malloc(n1 * n2 * sizeof(y[0][0]));
	if (!y) err(1, "Failed to allocate y for %d complex doubles", n1 * n2);

	complex double *w2 = compute_omegas(n2 / 2, flag);
	complex double *w1 = (n1 == n2) ? w2 : compute_omegas(n1 / 2, flag);


	/* FFT in y */
	for (int i = 0; i < n1; i ++)
		stockham_1D(n2, x[i], y[i], w2);

	/* FFT in x */
	cooley_tukey(n1, n2, x, w1);


	free(y);
	free(w2);
	if (n2 != n1) free(w1);
}

/* 3D Fourier transform:
   The index for x is mapped to (i,j,k) by x[i][j][k] = x[i*n2*n3 + j*n3 + k],
   i.e. the row major convention of C.  All indices start from 0.
   This algorithm requires working space of n2*n3.
   Stockham is efficient, good stride feature, but takes extra
   memory same size as input data; Cooley-Tukey is in place,
   so we take a compromise of the two.
*/
void seq_fft3D(int n1, int n2, int n3, complex double (*x)[n2][n3], int flag)
{
	complex double (*y)[n2][n3] = (complex double (*)[n2][n3])malloc(n1 * n2 * n3 * sizeof(y[0][0][0]));
	if (!y) err(1, "Failed to allocate y for %d complex doubles", n1 * n2 * n3);

	complex double *w3 = compute_omegas(n3 / 2, flag);
	complex double *w2 = (n2 == n3) ? w3 : compute_omegas(n2 / 2, flag);
	complex double *w1 = (n1 == n3 || n1 == n2) ? (n1 == n3 ? w3 : w2) : compute_omegas(n1 / 2, flag);


	/* FFT in z */
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			stockham_1D(n3, x[i][j], y[i][j], w3);

	/* FFT in y */
	for (int i = 0; i < n1; i++)
		stockham(n2, n3, x[i], y[i], w2);

	/* FFT in x */
	int n23 = n2 * n3;
	cooley_tukey(n1, n23, (complex double (*)[n23])x[0], w1);


	free(y);
	free(w3);
	if (n3 != n2) free(w2);
	if (n3 != n1 && n2 != n1) free(w1);
}
