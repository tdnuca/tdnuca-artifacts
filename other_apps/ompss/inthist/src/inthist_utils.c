#include "inthist_utils.h"


#include <math.h>
#include <stdio.h>
#include <string.h>

#include <unistd.h>
#include <stdlib.h>
#ifndef _PAGE_SIZE
#define _PAGE_SIZE sysconf(_SC_PAGESIZE)
#endif 

int util_alloc(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth],\
  unsigned int *vhalos[imheight],unsigned int *hhalos[imwidth],int bheight,int bwidth,int oblocksize,int binsize)
{
  int num_sockets;
  count_sockets(&num_sockets);
  
  int i;
  for(i=0;i<imheight;i++)
  {
    current_socket(i%num_sockets);
    int j;
    for(j=0;j<imwidth;j++)
    {

      if (img != NULL) {
        if(posix_memalign((void**) &img[i][j], _PAGE_SIZE, bheight*bwidth*sizeof(float))) {
	  fprintf(stderr,"out of memory!\n");
          return 1;
        }
      }

      if(posix_memalign((void**) &inthist[i][j], _PAGE_SIZE, oblocksize*sizeof(unsigned int))) {
	fprintf(stderr,"out of memory!\n");
	return 1;
      }

      if (img != NULL) {
#if USE_OMPSS
#pragma omp task out([bheight*bwidth](img[i][j])) firstprivate(i, j, bheight, bwidth) label(zero_img)
#endif
        memset(img[i][j], 0, bheight*bwidth*sizeof(float));
      }

#if USE_OMPSS
#pragma omp task out([oblocksize](inthist[i][j])) firstprivate(i, j, oblocksize) label(zero_hist)
#endif
      memset(inthist[i][j], 0, oblocksize*sizeof(unsigned int));
    }
  }

  for(i=0;i<imheight;i++)
  {
    current_socket(i%num_sockets);
    if(posix_memalign((void**) &vhalos[i], _PAGE_SIZE, binsize*bheight*sizeof(unsigned int)))
    {
      fprintf(stderr,"out of memory!\n");
      return 1;
    }

#if USE_OMPSS
#pragma omp task out([binsize*bheight](vhalos[i])) firstprivate(i) label(zero_vhalo)
#endif
    memset(vhalos[i], 0, binsize*bheight*sizeof(unsigned int));
  }

  for(i=0;i<imwidth;i++)
  {
    current_socket(i%num_sockets);
    if(posix_memalign((void**) &hhalos[i], _PAGE_SIZE, binsize*bwidth*sizeof(unsigned int)))
    {
      fprintf(stderr,"out of memory!\n");
      return 1;
    }

#if USE_OMPSS
#pragma omp task out([binsize*bwidth](hhalos[i])) firstprivate(i) label(zero_hhalo)
#endif
    memset(hhalos[i], 0, binsize*bwidth*sizeof(unsigned int));
  }

  
  return 0;
}



void util_freehalos(int imheight,int imwidth,unsigned int *vhalos[imheight],unsigned int *hhalos[imwidth])
{
  int i;
  for(i=0;i<imheight;i++)
  {
    free(vhalos[i]);
  }

  for(i=0;i<imwidth;i++)
  {
    free(hhalos[i]);
  }
}


void util_freeimage(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth])
{
  int maxdim=imwidth>imheight?imwidth:imheight;

  int i;
  for(i=0;i<imheight;i++)
  {
    int j;
    for(j=0;j<imwidth;j++)
    {
      if (img != NULL)
        free(img[i][j]);
      free(inthist[i][j]);
    }
  }
}


void util_hist(const char fname[],unsigned int *endhist,int histsize,float binstride)
{
  printf("creating image histogram...\n");
  FILE *f=fopen(fname,"w");

  float binstart;
  int i;
  for(i=0,binstart=0;i<histsize;i++,binstart+=binstride)
  {
    fprintf(f,"%f %i\n",binstart,endhist[i]);
  }
}


void util_reset(int imheight,int imwidth,unsigned int *vhalos[imheight],unsigned int *hhalos[imwidth],int bheight,int bwidth,int binsize)
{
  int num_sockets;
  count_sockets(&num_sockets);

  int i;
  for(i=0;i<imheight;i++)
  {
    current_socket(i%num_sockets);
#if USE_OMPSS
#pragma omp task out([binsize*bheight](vhalos[i]))
#endif
    memset(vhalos[i],0,binsize*bheight*sizeof(unsigned int));
  }

  for(i=0;i<imwidth;i++)
  {
    current_socket(i%num_sockets);
#if USE_OMPSS
#pragma omp task out([binsize*bwidth](vhalos[i]))
#endif
    memset(hhalos[i],0,binsize*bwidth*sizeof(unsigned int));
  }
}
