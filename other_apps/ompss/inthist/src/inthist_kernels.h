#ifndef __INT_HISTOGRAM_KERNELS__
#define __INT_HISTOGRAM_KERNELS__

#include "inthist_extra.h"

#define PIXEL_T float
#define BIN_T	unsigned int


#if USE_OMPSS
#pragma omp task in([bheight*bwidth]im)	\
  inout(whalo[0;binsize*bheight])	\
  out([oblocksize]inthist)		\
  label(hscan)
#endif
void inthist_hscan(int bheight,int bwidth,int oblocksize,int binsize,float binstride,\
  PIXEL_T *im,\
  BIN_T *whalo,\
  BIN_T *inthist);


#if USE_OMPSS
#pragma omp task inout([binsize*bwidth]nhalo) \
  inout([oblocksize]inthist)		      \
  label(vscan)
#endif
void inthist_vscan(int bheight,int bwidth,int oblocksize,int binsize,float binstride,\
  BIN_T *nhalo,\
  BIN_T *inthist);


#endif // __INT_HISTOGRAM_KERNELS__
