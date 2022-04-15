#include "inthist_kernels.h"

#include <stdio.h>
#include <math.h>

void warmup(int size,unsigned int *buffer,char *dummy)
{
}


void inthist_hscan(int bheight,int bwidth,int oblocksize,int binsize,float binstride,\
  PIXEL_T *im,
  BIN_T *whalo,\
  BIN_T *inthist)
{
  unsigned int *wstart=whalo;
  unsigned int *ehalo=whalo;
  unsigned int *w;

  int i;
  for(i=0;i<bheight;i++)
  {
    w=wstart;
    int j;
    for(j=0;j<bwidth;j++)
    {
      unsigned int *outstart=inthist;

      unsigned int border=(j==(bwidth-1));
      unsigned int incr=border?1:0;

      float pixel=im[i*bwidth+j];
      float bin_pr=pixel/binstride;
      float floored_bin_pr=floor(bin_pr);
      unsigned int floorbin=floored_bin_pr;

      int jj;
      for(jj=0;jj<binsize;jj++)
      {
	unsigned int binfit=floorbin==jj;
	*inthist=*w+binfit;

	*ehalo=border?*inthist:*ehalo;

	w++;
	inthist++;
	ehalo+=incr;
      }

      w=outstart;
    }
    wstart+=binsize;
  }
}


void inthist_vscan(int bheight,int bwidth,int oblocksize,int binsize,float binstride,\
  BIN_T *nhalo,\
  BIN_T *inthist)
{
  unsigned int *shalo=nhalo;

  int i;
  for(i=0;i<bheight-1;i++)
  {
    unsigned int *outstart=inthist;

    int j;
    for(j=0;j<bwidth;j++)
    {
      int jj;
      for(jj=0;jj<binsize;jj++)
      {
	*inthist=*nhalo+*inthist;
	
	nhalo++;
	inthist++;
      }
    }

    nhalo=outstart;
  }

  int j;
  for(j=0;j<bwidth;j++)
  {
    int jj;
    for(jj=0;jj<binsize;jj++)
    {
      *inthist=*nhalo+*inthist;
      *shalo=*inthist;
	
      shalo++;
      nhalo++;
      inthist++;
    }
  }
}

