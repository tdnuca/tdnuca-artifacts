#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#include "inthist_extra.h"
#include "inthist_check.h"
#include "inthist_image.h"
#include "inthist_kernels.h"
#include "inthist_utils.h"


void writestats(const char statfname[],double imcnt,double elapsed);

int main(int argc,char *argv[])
{
    if(argc<8)
    {
        fprintf(stderr,"use: %s imw imh bw bh bins binstride ims [outfile]\n",argv[0]);
        return 1;
    }

    int trueimwidth=atoi(argv[1]);
    int trueimheight=atoi(argv[2]);
    int bwidth=atoi(argv[3]);
    int bheight=atoi(argv[4]);
    int imwidth=(trueimwidth+bwidth-1)/bwidth;
    int imheight=(trueimheight+bheight-1)/bheight;
    int truebinsize=atoi(argv[5]);
    unsigned int binsize=(truebinsize+3)&~3;
    float binstride=strtof(argv[6],NULL);
    int imcnt=atoi(argv[7]);


    printf("image %ix%i blocks %ix%i binsize %i binstride %f imcnt %i\n",imwidth,imheight,bwidth,bheight,binsize,binstride,imcnt);

    int lasthistidx=(bheight-1)*bwidth*binsize+(bwidth-1)*binsize;

    // Set up the images and halos.
    float *im[imheight][imwidth];
#if SEQUENCE_MODE
    unsigned int *inthist[2][imheight][imwidth];
    unsigned int *vhalos[2][imheight];
    unsigned int *hhalos[2][imwidth];
#else
    unsigned int *inthist[imheight][imwidth];
    unsigned int *vhalos[imheight];
    unsigned int *hhalos[imwidth];
#endif

    int oblocksize=binsize*bwidth*bheight;
    // unsigned int bsize=bwidth*bheight;

#if SEQUENCE_MODE
    if(util_alloc(imheight,imwidth,im,inthist[0],vhalos[0],hhalos[0],bheight,bwidth,oblocksize,binsize))
    {
        fprintf(stderr,"Exiting...\n");
    }
    if(util_alloc(imheight,imwidth,NULL,inthist[1],vhalos[1],hhalos[1],bheight,bwidth,oblocksize,binsize))
    {
        fprintf(stderr,"Exiting...\n");
    }
#else
    if(util_alloc(imheight,imwidth,im,inthist,vhalos,hhalos,bheight,bwidth,oblocksize,binsize))
    {
        fprintf(stderr,"Exiting...\n");
    }
#endif

    int errorc=0;
    int imi;

    unsigned long elapsed=0;
    unsigned long fullelapsed=0;
    struct timeval start,stop;

#if SEQUENCE_MODE
    gettimeofday(&start,NULL);
#endif


    int num_sockets;
    count_sockets(&num_sockets);

    for(imi=0;imi<imcnt;imi++)
    {
        float max=readimg(imheight,imwidth,im,bheight,bwidth);

#if !SEQUENCE_MODE
        gettimeofday(&start,NULL);
#else
        int imflip=imi&1;
#endif

        int i;
        for(i=0;i<imheight;i++)
        {
            int j;
            current_socket(i%num_sockets);
            for(j=0;j<imwidth;j++)
            {
#if SEQUENCE_MODE
                inthist_hscan(bheight,bwidth,oblocksize,binsize,binstride,\
                        im[i][j],\
                        vhalos[imflip][i],\
                        inthist[imflip][i][j]);
#else 
                inthist_hscan(bheight,bwidth,oblocksize,binsize,binstride,\
                        im[i][j],\
                        vhalos[i],\
                        inthist[i][j]);
#endif
            }
        }

        for(i=0;i<imwidth;i++)
        {
            int j;
            for(j=0;j<imheight;j++)
            {
                current_socket(j%num_sockets);
#if SEQUENCE_MODE 
                inthist_vscan(bheight,bwidth,oblocksize,binsize,binstride,\
                        hhalos[imflip][i],\
                        inthist[imflip][j][i]);
#else
                inthist_vscan(bheight,bwidth,oblocksize,binsize,binstride,\
                        hhalos[i],\
                        inthist[j][i]);
#endif
            }
        }




#if !SEQUENCE_MODE
#if USE_OMPSS
#pragma omp taskwait 
#endif

        gettimeofday(&stop,NULL);
        unsigned long fullpart = 1000000 * (stop.tv_sec - start.tv_sec);
        fullpart += stop.tv_usec - start.tv_usec;
        elapsed+=fullpart;

#if INTHIST_CHECK
#if USE_OMPSS
#pragma omp taskwait 
#endif
        errorc+=inthist_check(imheight,imwidth,im,inthist,bheight,bwidth,binsize,binstride);
#endif
#endif
    }




#if SEQUENCE_MODE
#if USE_OMPSS // not needed, SEQUENCE_MODE implies OmpSs, but...
#pragma omp taskwait 
#endif
    gettimeofday(&stop,NULL);
    fullelapsed = 1000000 * (stop.tv_sec - start.tv_sec);
    fullelapsed += stop.tv_usec - start.tv_usec;
    elapsed = fullelapsed;
#endif

    double imcntfl=imcnt;
    double elapsedfl=elapsed;

    printf("inthist wall clock            : %lu us\n",elapsed);
    printf("inthist wall clock per image  : %.2f us\n",elapsedfl/imcntfl);
    printf("inthist frate                 : %.2f fps\n",imcntfl/(elapsedfl/1000000.));
    printf("inthist blockrate             : %.2f blps\n",(imcntfl*(imheight*imwidth))/(elapsedfl/1000000.));


    if (argc > 8)
        writestats(argv[8],imcntfl,elapsedfl);


#if SEQUENCE_MODE
    util_freehalos(imheight,imwidth,vhalos[0],hhalos[0]);
    util_freehalos(imheight,imwidth,vhalos[1],hhalos[1]);
    util_freeimage(imheight,imwidth,im,inthist[0]);
    util_freeimage(imheight,imwidth,NULL,inthist[1]);
#else
    util_freehalos(imheight,imwidth,vhalos,hhalos);
#endif

    fprintf(stderr, "ERROR = %d\n", errorc);

    return errorc;
}



void writestats(const char statfname[],double imcnt,double elapsed)
{
    FILE *statfile;
    if(statfile=fopen(statfname,"w"))
    {
        fprintf(statfile,"inthist crossweave frate : %.2f\n",imcnt/(elapsed/1000000.));
        fprintf(statfile,"inthist crossweave image : %.2f\n",elapsed/imcnt);
        fprintf(statfile,"inthist crossweave full  : %.2f\n",elapsed);
    }
    fclose(statfile);
}

