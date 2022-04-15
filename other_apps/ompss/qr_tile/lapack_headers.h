#ifndef __LAPACK_HEADERS_H
#define __LAPACK_HEADERS_H

#ifdef __cplusplus
extern "C" {
#endif

  typedef int Integer;

  int blocksize(int bs);

  int ilaenv_(int, char *, char *, int, int, int, int);

  void dgeqrt_(Integer *m, Integer *n, Integer *nb,
	       double *A, Integer *lda,
	       double *T, Integer *ldt,
	       double *work,
	       Integer *info);

  void dlarfb_(const char *side, const char *trans, const char *direct, const char *storev,
	       Integer *m, Integer *n, Integer *k,
	       double *V, Integer *ldv,
	       double *T, Integer *ldt,
	       double *C, Integer *ldc,
	       double *work, Integer *ldwork);

  void dtpqrt_(Integer *m, Integer *n, Integer *l, Integer *nb,
	       double *A, Integer *lda,
	       double *B, Integer *ldb,
	       double *T, Integer *ldt,
	       double *work,
	       Integer *info);

  void dtpmqrt_(const char *side, const char *trans,
		Integer *m, Integer *n, Integer *k, Integer *l, Integer *nb,
		double *V, Integer *ldv,
		double *T, Integer *ldt,
		double *A, Integer *lda,
		double *B, Integer *ldb,
		double *work,
		Integer *info);

#ifdef __cplusplus
}
#endif

#endif
