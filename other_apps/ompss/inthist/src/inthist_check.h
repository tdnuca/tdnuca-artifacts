#ifndef __INTHIST_CHECK_H__
#define __INTHIST_CHECK_H__

#include "inthist_extra.h"

#if INTHIST_CHECK
int inthist_check(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth],\
  int bheight,int bwidth,int binsize,float binstride);
#else
static inline __attribute__((always_inline)) int inthist_check(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth],\
  int bheight,int bwidth,int binsize,float binstride)
{
  return 0;
}
#endif


#endif // __INTHIST_CHECK_H__
