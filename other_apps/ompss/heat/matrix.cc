#include "matrix.hh"
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <unistd.h>

#ifndef _NANOS_H_
static void nanos_current_socket(int x) {}
static void nanos_get_num_sockets(int *x)
{
  *x = 1;
}
#endif

#ifndef INIT_TASKS
#define INIT_TASKS 1
#endif

#if INIT_TASKS
#pragma omp task label(border_step) inout([BS]b)
#endif
static void border_step(const long BS,
			const double limit,
			const double range,
			const double temp,
		        long x,
			const double src_x,
			const double y_sq,
			double *b);

#if INIT_TASKS
#pragma omp task label(fill_zero) inout([size]b)
#endif
static void fill_zero(const long size, double *b);

Matrix::Matrix(unsigned long NB, unsigned long BS, bool dbuffer)
  : NB(NB),
    BS(BS),
    socket_change(NB),
    dbuffer(dbuffer),
    current(0)
{
  nanos_get_num_sockets(&num_sockets);
  u[0] = new double*[NB*NB];

#if INIT_TASKS
  socket_change = (NB + num_sockets - 1)/num_sockets;
#endif
    
  for (unsigned long i = 0; i < NB; ++i) {
    for (unsigned long j = 0; j < NB; ++j) {
      int res = posix_memalign((void**) &u[0][i*NB+j],
			       sysconf(_SC_PAGESIZE),
			       BS*BS*sizeof(double));
      if (res != 0) {
	std::cerr << "There was an error allocating matrix `u'!" << std::endl;
      }

#if INIT_TASKS
      nanos_current_socket(j/socket_change);
#endif
      fill_zero(BS*BS, u[0][i*NB+j]);
    }
  }

  vertical[0] = new double*[NB*(NB + 1)];
  for (unsigned long i = 0; i < NB; ++i) {
    for (unsigned long j = 0; j < NB + 1; ++j) {

      int res = posix_memalign((void **) &vertical[0][i*(NB+1)+j],
			       sysconf(_SC_PAGESIZE),
			       BS*sizeof(double));
      if (res != 0) {
	std::cerr << "There was an error allocating vertical halo!" << std::endl;
      }
	
#if INIT_TASKS
      nanos_current_socket((j > i ? j-1 : j)/socket_change);
#endif
      fill_zero(BS, vertical[0][i*(NB+1)+j]);
    }
  }

  horizontal[0] = new double*[(NB + 1)*NB];
  for (unsigned long i = 0; i < NB + 1; ++i) {
    for (unsigned long j = 0; j < NB; ++j) {
      int res = posix_memalign((void **) &horizontal[0][i*NB+j],
			       sysconf(_SC_PAGESIZE),
			       BS*sizeof(double));
      if (res != 0) {
	std::cerr << "There was an error allocating horizontal halo!" << std::endl;
      }
	
#if INIT_TASKS
      nanos_current_socket(j/socket_change);
#endif
      fill_zero(BS, horizontal[0][i*NB+j]);
    }
  }

  if (dbuffer) {
    u[1] = new double*[NB*NB];
    for (unsigned long i = 0; i < NB; ++i) {
      for (unsigned long j = 0; j < NB; ++j) {
	int res = posix_memalign((void**) &u[1][i*NB+j],
				 sysconf(_SC_PAGESIZE),
				 BS*BS*sizeof(double));
	if (res != 0) {
	  std::cerr << "There was an error allocating matrix `u'!" << std::endl;
	}

// These inits are not needed because they will be overwritten before being used.
// Left them just for clarity.
// #if INIT_TASKS
// 	nanos_current_socket(j/socket_change);
// #endif
// 	fill_zero(BS*BS, u[1][i*NB+j]);
      }
    }

    vertical[1] = new double*[NB*(NB + 1)];
    for (unsigned long i = 0; i < NB; ++i) {
      vertical[1][i*(NB+1)+0] = vertical[0][i*(NB+1)+0];
      vertical[1][i*(NB+1)+NB] = vertical[0][i*(NB+1)+NB];
      for (unsigned long j = 1; j < NB; ++j) {
	int res = posix_memalign((void **) &vertical[1][i*(NB+1)+j],
				 sysconf(_SC_PAGESIZE),
				 BS*sizeof(double));
	if (res != 0) {
	  std::cerr << "There was an error allocating vertical halo!" << std::endl;
	}
	
#if INIT_TASKS
	nanos_current_socket((j > i ? j-1 : j)/socket_change);
#endif
	fill_zero(BS, vertical[1][i*(NB+1)+j]);
      }
    }

    horizontal[1] = new double*[(NB + 1)*NB];
    for (unsigned long j = 0; j < NB; ++j) {
      horizontal[1][0*NB+j] = horizontal[0][0*NB+j];
      horizontal[1][NB*NB+j] = horizontal[0][NB*NB+j];
      for (unsigned long i = 1; i < NB; ++i) {
	int res = posix_memalign((void **) &horizontal[1][i*NB+j],
				 sysconf(_SC_PAGESIZE),
				 BS*sizeof(double));
	if (res != 0) {
	  std::cerr << "There was an error allocating horizontal halo!" << std::endl;
	}
	
#if INIT_TASKS
	nanos_current_socket(j/socket_change);
#endif
	fill_zero(BS, horizontal[1][i*NB+j]);
      }
    }
  }
}

Matrix::~Matrix()
{
  for (unsigned long i = 0; i < NB*NB; ++i) {
    free(u[0][i]);
  }
  delete[] u[0];
    
  for (unsigned long i = 0; i < NB; ++i) {
    for (unsigned long j = 0; j < NB + 1; ++j) {
      free(vertical[0][i*(NB+1)+j]);
    }
  }
  delete[] vertical[0];

  for (unsigned long i = 0; i < NB + 1; ++i) {
    for (unsigned long j = 0; j < NB; ++j) {
      free(horizontal[0][i*NB+j]);
    }
  }
  delete[] horizontal[0];

  if (dbuffer) {
    for (unsigned long i = 0; i < NB*NB; ++i) {
      free(u[1][i]);
    }
    delete[] u[1];
    
    for (unsigned long i = 0; i < NB; ++i) {
      for (unsigned long j = 1; j < NB; ++j) {
	free(vertical[1][i*(NB+1)+j]);
      }
    }
    delete[] vertical[1];

    for (unsigned long i = 1; i < NB; ++i) {
      for (unsigned long j = 0; j < NB; ++j) {
	free(horizontal[1][i*NB+j]);
      }
    }
    delete[] horizontal[1];
  }
}

void Matrix::setBorder(const heatsrc_t &src)
{
  double limit = NB*BS-1;

  // TOP row
  double y_sq = std::pow(src.posy, 2.0);
  for (unsigned long bl = 0; bl < NB; ++bl) {
#if INIT_TASKS
    nanos_current_socket(bl/socket_change);
#endif
    border_step(BS, limit, src.range, src.temp,
		bl*BS, src.posx, y_sq,
		top(0, bl));
  }
    
  // RIGHT column
  y_sq = std::pow(1.0 - src.posx, 2.0);
  for (unsigned long bl = 0; bl < NB; ++bl) {
#if INIT_TASKS
    nanos_current_socket(num_sockets - 1);
#endif
    border_step(BS, limit, src.range, src.temp,
		bl*BS, src.posy, y_sq,
		right(bl, NB-1));
  }
    
  // BOTTOM row
  y_sq = std::pow(1.0 - src.posy, 2.0);
  for (unsigned long bl = 0; bl < NB; ++bl) {
#if INIT_TASKS
    nanos_current_socket(bl/socket_change);
#endif
    border_step(BS, limit, src.range, src.temp,
		bl*BS, src.posx, y_sq,
		bottom(NB-1, bl));
  }
    
  // LEFT column
  y_sq = std::pow(src.posx, 2.0);
  for (unsigned long bl = 0; bl < NB; ++bl) {
#if INIT_TASKS
    nanos_current_socket(0);
#endif
    border_step(BS, limit, src.range, src.temp,
		bl*BS, src.posy, y_sq,
		left(bl, 0));
  }
}


static void fill_zero(const long size, double *b)
{
  std::fill(&b[0], &b[size], 0.0);
}


static void border_step(const long BS,
			const double limit,
			const double range,
			const double temp,
		        long x,
			const double src_x,
			const double y_sq,
			double *b)
{
  for (long i = 0; i < BS; ++i) {
    double dist = std::sqrt(std::pow(double(x+i)/limit - src_x, 2.0) + y_sq);
    
    if (dist <= range)
      b[i] += (range-dist) / range * temp;
  }
}
