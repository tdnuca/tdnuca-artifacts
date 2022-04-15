#ifndef __INT_HISTOGRAM_UTILS_H__
#define __INT_HISTOGRAM_UTILS_H__

#include "inthist_extra.h"

int util_alloc(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth],\
  unsigned int *vhalos[imheight],unsigned int *hhalos[imwidth],int bheight,int bwidth,int oblocksize,int binsize);

void util_freehalos(int imheight,int imwidth,unsigned int *vhalos[imheight],unsigned int *hhalos[imwidth]);

void util_freeimage(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth]);

void util_hist(const char fname[],unsigned int *endhist,int histsize,float binstride);

void util_reset(int imheight,int imwidth,unsigned int *vhalos[imheight],unsigned int *hhalos[imwidth],int bheight,int bwidth,int binsize);

#endif // __INT_HISTOGRAM_UTILS_H__
