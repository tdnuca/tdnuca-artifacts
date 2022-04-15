#include "matrix.hh"
#include "heat.hh"
#include "algorithms_inner.hh"

#ifndef _NANOS_H_
static void nanos_current_socket(int x) {}
static void nanos_get_num_sockets(int *x)
{
  *x = 1;
}
#endif

#pragma omp task label(heat_step)		\
  inout( [BS*BS]mu )				\
  inout( [BS]top,				\
	 [BS]bottom,				\
	 [BS]left,				\
	 [BS]right )				\
  shared(*sum)
static void heat_step(const unsigned long BS,
		      double *mu,
		      double *top, const bool wr_top,
		      double *bottom, const bool wr_bottom,
		      double *left, const bool wr_left,
		      double *right, const bool wr_right,
		      double *sum, const bool check);

#pragma omp task label(heat_step)		\
  in( [BS*BS]oldmu )				\
  out( [BS*BS]mu )				\
  inout( [BS]top,				\
	 [BS]bottom,				\
	 [BS]left,				\
	 [BS]right )				\
  shared(*sum)
static void heat_step(const unsigned long BS,
                      double *oldmu,
		      double *mu,
		      double *top, const bool wr_top,
		      double *bottom, const bool wr_bottom,
		      double *left, const bool wr_left,
		      double *right, const bool wr_right,
		      double *sum, const bool check);


void relax_gauss(Matrix &m, bool check, double *residual)
{
  if (check)
    *residual = 0.0;
  
  long BS = m.BS;
  long NB = m.NB;
  long max_diag = 2*NB - 1;
  for (long diag = 0; diag < max_diag; ++diag) {
    long start_i = std::max(0L, diag - NB + 1);
    long max_i = std::min(diag + 1, NB);
    for (long i = start_i; i < max_i; ++i) {
      long j = diag - i;
      nanos_current_socket(j/m.socket_change);
      heat_step(BS,
		m.block(i, j),
		m.top(i, j), i != 0,
		m.bottom(i, j), i != NB - 1,
		m.left(i, j), j != 0,
		m.right(i, j), j != NB - 1,
		residual, check);
    }
  }
}

void relax_jacobi(Matrix &m, bool check, double *residual)
{
  if (check)
    *residual = 0.0;

  ++m; // change of matrix!
  long BS = m.BS;
  long NB = m.NB;
  for (long i = 0; i < NB; ++i) {
    for (long j = 0; j < NB; ++j) {
      nanos_current_socket(j/m.socket_change);
      heat_step(BS,
		m.oldblock(i, j),
		m.block(i, j),
		m.top(i, j), i != 0,
		m.bottom(i, j), i != NB - 1,
		m.left(i, j), j != 0,
		m.right(i, j), j != NB - 1,
		residual, check);
    }
  }
}


void relax_redblack(Matrix &m, bool check, double *residual)
{
  if (check)
    *residual = 0.0;
  
  long BS = m.BS;
  long NB = m.NB;
  for (long rb = 0; rb < 2; ++rb) {
    for (long j = 0; j < NB; ++j) {
      long i_start = (j + rb)%2;
      for (long i = i_start; i < NB; i += 2) {
        nanos_current_socket(j/m.socket_change);
	heat_step(BS,
		  m.block(i, j),
		  m.top(i, j), i != 0,
		  m.bottom(i, j), i != NB - 1,
		  m.left(i, j), j != 0,
		  m.right(i, j), j != NB - 1,
		  residual, check);
      }
    }
  }
}


static void heat_step(const unsigned long BS,
		      double *mu,
		      double *top, const bool wr_top,
		      double *bottom, const bool wr_bottom,
		      double *left, const bool wr_left,
		      double *right, const bool wr_right,
		      double *sum, const bool check)
{
  double my_sum = inner_heat_step(BS,
				  mu,
				  top, wr_top,
				  bottom, wr_bottom,
				  left, wr_left,
				  right, wr_right,
				  check);
  if (check) {
#pragma omp atomic
    *sum += my_sum;
  }
}

static void heat_step(const unsigned long BS,
		      double *oldmu,
		      double *mu,
		      double *top, const bool wr_top,
		      double *bottom, const bool wr_bottom,
		      double *left, const bool wr_left,
		      double *right, const bool wr_right,
		      double *sum, const bool check)
{
  double my_sum = inner_heat_step(BS,
				  oldmu,
				  mu,
				  top, wr_top,
				  bottom, wr_bottom,
				  left, wr_left,
				  right, wr_right,
				  check);
  if (check) {
#pragma omp atomic
    *sum += my_sum;
  }
}

