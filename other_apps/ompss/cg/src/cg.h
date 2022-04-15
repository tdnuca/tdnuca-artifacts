#ifndef CG_H_INCLUDED
#define CG_H_INCLUDED

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "global.h"
#include "matrix.h"

void solve_cg(const Matrix *A, const double *b, double *iterate, double convergence_thres);

#ifdef USE_MPI
void determine_mpi_neighbours(const Matrix *A, const int from_row, const int to_row, const int mpi_size, int *first, int *last);
void setup_exchange_vect(const int mpi_rank, const int first, const int last, const int tag, double *v, MPI_Request v_req[]);
#endif

// all the algorithmical steps of CG that will be subdivided into tasks :
void update_gradient(double *gradient, double *Ap, double *alpha);
void recompute_gradient_mvm(const Matrix *A, double *iterate, double *Aiterate);
void recompute_gradient_update(double *gradient, double *Aiterate, const double *b);
void update_p(double *p, double *old_p, double *gradient, double *beta);
void update_iterate(double *iterate, double *p, double *alpha);
void compute_Ap(const Matrix *A, double *p, double *Ap);

void scalar_product_task(const double *p, const double *Ap, double* r);
void norm_task(const double *v, double* r);

#ifndef USE_MPI
void compute_beta(const double *err_sq, const double *old_err_sq, double *beta);
#else
void compute_beta(double *err_sq, const double *old_err_sq, double *beta);
#endif
void compute_alpha(double *err_sq, double *normA_p_sq, double *old_err_sq, double *alpha);

#endif // CG_H_INCLUDED
