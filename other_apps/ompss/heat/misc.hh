#ifndef _MISC_HH
#define _MISC_HH
#include "matrix.hh"
#include <cstdio>
#include <iostream>

struct algoparam_t
{
  unsigned maxiter;       // maximum number of iterations
  unsigned resolution;    // spatial resolution
  int algorithm;          // 0=>Jacobi, 1=>Gauss
  Matrix *m;
  unsigned visres;        // visualization resolution
  double *uvis;
  unsigned numsrcs;     // number of heat sources
  heatsrc_t *heatsrcs;

  algoparam_t();
  ~algoparam_t();
  int read_input( std::istream& infile );
  int init(unsigned long NB);
  void print_params() const;
  void write_image( std::FILE * f );
  int coarsen();
};

double wtime();

struct color_t
{
  int r, g, b;
};

typedef color_t colormap_t[256];

extern colormap_t magma; 

#endif
