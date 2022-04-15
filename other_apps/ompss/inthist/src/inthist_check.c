#include "inthist_check.h"


#if INTHIST_CHECK


#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int inthist_check(int imheight,int imwidth,float *img[imheight][imwidth],unsigned int *inthist[imheight][imwidth],\
  int bheight,int bwidth,int binsize,float binstride)
{
  if(INTHIST_CHECK>1)
  {
    printf("checking...");
    fflush(0);
  }

  int maxdim=imwidth>imheight?imwidth:imheight;
  unsigned int *zerobin=(unsigned int*)calloc(binsize,sizeof(unsigned int));

  unsigned int errorc=0;

  unsigned int *prevnhalo=NULL;;

  int bi;
  for(bi=0;bi<imheight;bi++)
  {
    int topblock=bi==0;

    int ei;
    for(ei=0;ei<bheight;ei++)
    {
      int *whalo=zerobin;

      int toprow=ei==0;

      int bj;
      for(bj=0;bj<imwidth;bj++)
      {
	unsigned int *inthistblock=inthist[bi][bj];
	inthistblock+=bwidth*binsize*ei;

	float *imblock=img[bi][bj];
	imblock+=bwidth*ei;

	int leftblock=bj==0;

	int ej;
	for(ej=0;ej<bwidth;ej++)
	{
	  int leftcol=ej==0;

	  int *nextwhalo=inthistblock;


	  float pixel=imblock[ej];

	  float bin_pr=pixel/binstride;
	  float floored_bin_pr=floor(bin_pr);
	  unsigned int floorbin=floored_bin_pr;


	  int *nblock=inthist[bi-1][bj];
	  nblock+=bwidth*binsize*(bheight-1);
	  nblock+=ej*binsize;
	  
	  unsigned int *nhalo=inthistblock-binsize*bwidth;
	  nhalo=(topblock && toprow)?zerobin:nhalo;
	  nhalo=(!topblock && toprow)?nblock:nhalo;


	  int *nwhalo=(leftblock && leftcol)?zerobin:prevnhalo;

	  prevnhalo=nhalo;
	  
	  int bini;
	  for(bini=0;bini<binsize;bini++)
	  {
	    unsigned int binfit=floorbin==bini;
	    unsigned int binexp=binfit + *nhalo + *whalo - *nwhalo;

	    int binval=*inthistblock;

	    if(binval!=binexp)
	    {
	      errorc++;
	      if(INTHIST_CHECK>2)
	      {
		printf("\nb(%i,%i) e(%i,%i)=%f bin %i: exp %i found %i\n",bi,bj,ei,ej,pixel,bini,binexp,binval);
		printf("\t%i=%i+%i+%i-%i\n",binexp,binfit,*nhalo,*whalo,*nwhalo);
	      }
	    }

	    whalo++;
	    nhalo++;
	    nwhalo++;

	    inthistblock++;
	  }

	  whalo=nextwhalo;
	  nextwhalo=inthistblock;
	}
      }
    }
  }


  if(INTHIST_CHECK>1)
  {
    printf("done\n");
  }


  if(errorc)
  {
    printf("found %i errors in the histogram\n",errorc);
    return 1;
  }


  return 0;
}


#endif
