#include "inthist_image.h"
#include <stdlib.h>

float readimg(int imheight,int imwidth,float *img[imheight][imwidth],int bheight,int bwidth)
{
    int num_sockets;
    count_sockets(&num_sockets);

    int i;
    for(i=0;i<imheight;i++) {
        current_socket(i%num_sockets);
        int j;
        for(j=0;j<imwidth;j++) {
            float *bl=img[i][j];

#if USE_OMPSS
#pragma omp task inout([bheight*bwidth](img[i][j])) firstprivate(i, j, bheight, bwidth) label(read_rnd_block)
#endif
            {
                struct drand48_data myseed;
#if USE_OMPSS
                nanos_wd_t wd = nanos_current_wd();
                srand48_r(nanos_get_wd_id(wd), &myseed);
#else
                srand48_r(rand(), &myseed);
#endif
                int k;
                for(k=0;k<bheight;k++) {
                    int l;
                    for(l=0;l<bwidth;l++) {
                        // Random between 0 and 255
                        long int result;
                        lrand48_r(&myseed, &result);
                        bl[k*bwidth+l] = ((int) result) & 0xff;
                    }
                }
            }
        }
    }

    return 1;
}
