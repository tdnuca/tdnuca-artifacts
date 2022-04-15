#ifndef STENCIL_KERNEL
#define STENCIL_KERNEL

#ifndef RADIUS
#define RADIUS 2
#endif

void stencil_kernel_star( int n, double (* restrict in)[n][n], double (* restrict out)[n][n], int jt, int it, int tile_size, double (* restrict weight)[2*RADIUS+1][2*RADIUS+1] );

void stencil_kernel_star_3in( int n, double (* restrict in_center)[n][n],
                                     double (* restrict in_up)[n][n],
                                     double (* restrict in_down)[n][n],
                                     double (* restrict out)[n][n],
                                     int jt, int it, int tile_size, double (* restrict weight)[2*RADIUS+1][2*RADIUS+1] );
void stencil_kernel_star_3in_dump_in( int n, double (* restrict in_center)[n][n],
                                     double (* restrict in_up)[n][n],
                                     double (* restrict in_down)[n][n],
                                     int jt, int it, int tile_size, double (* restrict weight)[2*RADIUS+1][2*RADIUS+1] );
#endif
