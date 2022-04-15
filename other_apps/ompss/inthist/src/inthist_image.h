#ifndef __INTHIST_IMAGE_H__
#define __INTHIST_IMAGE_H__


#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "inthist_extra.h"

#if USE_RANDOMIMAGE

float readimg(int imheight,int imwidth,float *img[imheight][imwidth],int bheight,int bwidth);

#else

static inline __attribute__((always_inline)) float readimg(int imheight,int imwidth,float *img[imheight][imwidth],int bheight,int bwidth) 
{return 1;}

#endif


#endif // __INTHIST_IMAGE_H__
