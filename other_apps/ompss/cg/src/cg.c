#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "global.h"
#include "debug.h"
#include "counters.h"

#include "cg.h"

#if VERBOSE >= SHOW_TASKINFO
void *first_p = NULL;
static inline int which_p_copy(const void *ptr) { return (ptr == first_p ? 0 : 1); };
#endif

#ifdef USE_MPI
void determine_mpi_neighbours(const Matrix *A, const int from_row, const int to_row, const int mpi_size, int *first, int *last)
{
	// find our furthest neighbour in both direction
	int max_col = 0, min_col = A->n, i;
	for(i=from_row; i<to_row; i++)
	{
		// using the fact that each A->c is sorted on A->r[i]..A->r[i+1]-1
		if( A->c[A->r[i]] < min_col )
			min_col = A->c[A->r[i]];

		if( A->c[A->r[i+1]-1] > max_col )
			max_col = A->c[A->r[i+1]-1];
	}

	// now check with which MPI block we communicate, depending on its rows
	// this is reflexive since the matrix is symmetric
	for(i=0; i<mpi_size; i++)
	{
		if (min_col >= mpi_zonestart[i])
			continue;
		else if(min_col < mpi_zonestart[i] + mpi_zonesize[i])
			*first = i;
		if (max_col >= mpi_zonestart[i] && max_col < mpi_zonestart[i] + mpi_zonesize[i])
		{
			*last = i;
			break;
		}
	}

	#if VERBOSE >= SHOW_FAILINFO
	printf("mpi_rank %d spans columns %d to %d, thus communicating with ranks [%d..%d]\n", mpi_rank, min_col, max_col, *first, *last);
	#endif
}

void setup_exchange_vect(const int mpi_rank, const int first, const int last, const int tag, double *v, MPI_Request v_req[])
{
	int i, j=0;

	for(i=first; i<=last; i++)
		if(i!=mpi_rank)
			// recvs are always from i start,size
			MPI_Recv_init(v + mpi_zonestart[   i    ], mpi_zonesize[   i    ], MPI_DOUBLE, i, tag, MPI_COMM_WORLD, v_req+(j++));

	for(i=first; i<=last; i++)
		if(i!=mpi_rank)
			// sends are always from mpi_rank start,size
			MPI_Send_init(localrank_ptr(v), mpi_zonesize[mpi_rank], MPI_DOUBLE, i, tag, MPI_COMM_WORLD, v_req+(j++));
}
#endif

void scalar_product_task(const double *p, const double *Ap, double* r)
{
	int i;
	for(i=0; i < nb_blocks; i ++ )
	{
		int s = get_block_start(i), e = get_block_end(i);

		// r <- <p, Ap>
#ifndef USE_REDUCTIONS
		#pragma omp task concurrent(*r) in(p[s:e-1], Ap[s:e-1]) firstprivate(s, e) label(dotp) priority(10) no_copy_deps
		{
			double local_r = 0;
			int k;
			for(k=s; k < e; k++)
				local_r += p[k] * Ap[k];

			#pragma omp atomic
				*r += local_r;

			log_err(SHOW_TASKINFO, "Blockrow scalar product <p[%d], Ap> block %d finished = %e\n", which_p_copy(p), world_block(i), local_r);
		}
#else
		#pragma omp task reduction(+:[1]r) in(p[s:e-1], Ap[s:e-1]) firstprivate(s, e) label(dotp) priority(10) no_copy_deps
		{
			int k;
			for(k=s; k < e; k++)
				*r += p[k] * Ap[k];

			log_err(SHOW_TASKINFO, "Blockrow scalar product <p[%d], Ap> block %d finished = %e\n", which_p_copy(p), world_block(i), *r);
		}
#endif
	}
}

void norm_task(const double *v, double* r)
{
	int i;
	for(i=0; i < nb_blocks; i ++ )
	{
		int s = get_block_start(i), e = get_block_end(i);

		// r <- || v ||
#ifndef USE_REDUCTIONS
		#pragma omp task concurrent(*r) in(v[s:e-1]) firstprivate(s, e) label(norm) priority(10) no_copy_deps
		{
			double local_r = 0;
			int k;
			for(k=s; k<e; k++)
				local_r += v[k] * v[k];

			#pragma omp atomic
				*r += local_r;

			log_err(SHOW_TASKINFO, "Blockrow square norm || g || part %d finished = %e\n", world_block(i), local_r);
		}
#else
		#pragma omp task reduction(+:[1]r) in(v[s:e-1]) firstprivate(s, e) label(norm) priority(10) no_copy_deps
		{
			int k;
			for(k=s; k<e; k++)
				*r += v[k] * v[k];

			log_err(SHOW_TASKINFO, "Blockrow square norm || g || part %d finished = %e\n", world_block(i), *r);
		}
#endif
	}
}

void update_gradient(double *gradient, double *Ap, double *alpha)
{
	int i;
	for(i=0; i < nb_blocks; i++)
	{
		int s = get_block_start(i), e = get_block_end(i);

		#pragma omp task in(*alpha, Ap[s:e-1]) inout(gradient[s:e-1]) firstprivate(s, e) label(update_gradient) priority(10) no_copy_deps
		{
			int k;
			for(k=s; k<e; k++)
				gradient[k] -= (*alpha) * Ap[k];

			log_err(SHOW_TASKINFO, "Updating gradient part %d finished = %e with alpha = %e\n", world_block(i), norm(e-s, &(gradient[s])), *alpha);
		}
	}
}

void recompute_gradient_mvm(const Matrix *A, double *iterate, double *Aiterate)
{
	int i;
	for(i=0; i < nb_blocks; i++)
	{
		int s = get_block_start(i), e = get_block_end(i);

#ifndef USE_MPI
		// Aiterate <- A * iterate
		#pragma omp task in({iterate[get_block_start(b):get_block_end(b)-1], b=0:nb_blocks-1}) out(Aiterate[s:e-1]) firstprivate(s, e) label(AxIt) priority(10) no_copy_deps
#else
		// Aiterate <- A * iterate, here iterate is globally indexed
		#pragma omp task in({iterate[glob_block_start(blk):glob_block_end(blk)-1], blk=0:nb_glob_blocks-1}) out(Aiterate[s:e-1]) firstprivate(s, e) label(AxIt) priority(10) no_copy_deps
#endif
		{
			int k, l;
			for(l=s; l<e; l++)
			{
				Aiterate[l] = 0;

				for(k=A->r[l]; k < A->r[l+1]; k++)
					Aiterate[l] += A->v[k] * iterate[ A->c[k] ];
			}

			log_err(SHOW_TASKINFO, "A * x part %d finished = %e\n", world_block(i), norm(e-s, &(Aiterate[s])));
		}
	}
}

void recompute_gradient_update(double *gradient, double *Aiterate, const double *b)
{
	int i;
	for(i=0; i < nb_blocks; i++)
	{
		int s = get_block_start(i), e = get_block_end(i);

		// gradient <- b - Aiterate
		#pragma omp task in(b[s:e-1], Aiterate[s:e-1]) out(gradient[s:e-1]) firstprivate(s, e) label(b-AxIt) priority(10) no_copy_deps
		{
			int k;
			for(k=s; k<e; k++)
				gradient[k] = b[k] - Aiterate[k] ;

			log_err(SHOW_TASKINFO, "b - Ax part %d finished = %e\n", world_block(i), norm(e-s, &(gradient[s])));
		}
	}
}

void update_p(double *p, double *old_p, double *gradient, double *beta)
{
	int i;
	for(i=0; i < nb_blocks; i++)
	{
		int s = get_block_start(i), e = get_block_end(i);

		// p <- beta * old_p + gradient
		#pragma omp task in(*beta, gradient[s:e-1], old_p[s:e-1]) out(p[s:e-1]) firstprivate(s, e) label(update_p) priority(10) no_copy_deps
		{
			int k;
			for (k=s; k<e; k++)
				p[k] = (*beta) * old_p[k] + gradient[k];

			log_err(SHOW_TASKINFO, "Updating p[%d from %d] part %d finished = %e with beta = %e\n", which_p_copy(p), which_p_copy(old_p), world_block(i), norm(e-s, &(p[s])), *beta);
		}
	}
}

void compute_Ap(const Matrix *A, double *p, double *Ap)
{
	int i;
	for(i=0; i < nb_blocks; i++)
	{
		int s = get_block_start(i), e = get_block_end(i);

#ifndef USE_MPI
		// Ap <- A * p
		#pragma omp task in({p[get_block_start(b):get_block_end(b)-1], b=0:nb_blocks-1}) out(Ap[s:e-1]) firstprivate(s, e) label(Axp) priority(20) no_copy_deps
#else
		// Ap <- A * p, here p is globally indexed
		#pragma omp task in({p[glob_block_start(blk):glob_block_end(blk)-1], blk=0:nb_glob_blocks-1}) out(Ap[s:e-1]) firstprivate(s, e) label(Axp) priority(20) no_copy_deps
#endif
		{
			int k, l;
			for(l=s; l<e; l++)
			{
				Ap[l] = 0;

				for(k=A->r[l]; k < A->r[l+1]; k++)
					Ap[l] += A->v[k] * p[ A->c[k] ];
			}

			log_err(SHOW_TASKINFO, "A * p[%d] part %d finished = %e\n", which_p_copy(localrank_ptr(p)), world_block(i), norm(e-s, &(Ap[s])));
		}
	}
}

void update_iterate(double *iterate, double *p, double *alpha)
{
	int i;
	for(i=0; i < nb_blocks; i++)
	{
		int s = get_block_start(i), e = get_block_end(i);

		// iterate <- iterate - alpha * p
		#pragma omp task in(*alpha, p[s:e-1]) inout(iterate[s:e-1]) firstprivate(s, e) label(update_iterate) priority(5) no_copy_deps
		{
			int k;
			for(k=s; k<e; k++)
				iterate[k] += (*alpha) * p[k];

			log_err(SHOW_TASKINFO, "Updating it (from p[%d]) part %d finished = %e with alpha = %e\n", which_p_copy(p), world_block(i), norm(e-s, &(iterate[s])), *alpha);
		}
	}
}

// slight semantic difference between both versions here:
// since compute_alpha does the old_err <- err swap, err needs to be inout to account for the neighbour values
#ifndef USE_MPI
#pragma omp task in(*err_sq, *old_err_sq) out(*beta) label(compute_beta) priority(100) no_copy_deps
void compute_beta(const double *err_sq, const double *old_err_sq, double *beta)
#else
#pragma omp task in(*old_err_sq) out(*beta) inout(*err_sq) label(compute_beta) priority(100) no_copy_deps
void compute_beta(double *err_sq, const double *old_err_sq, double *beta)
#endif
{
#ifdef USE_MPI
	double loc_err_sq = *err_sq;
	MPI_Allreduce(&loc_err_sq, err_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

	// on first iterations of a (re)start, old_err_sq should be INFINITY so that beta = 0
	*beta = *err_sq / *old_err_sq;

	log_err(SHOW_TASKINFO, "Computing beta finished : err_sq = %e ; old_err_sq = %e ; beta = %e\n", *err_sq, *old_err_sq, *beta);
}

#pragma omp task inout(*normA_p_sq, *err_sq) out(*alpha, *old_err_sq) label(compute_alpha) priority(100) no_copy_deps
void compute_alpha(double *err_sq, double *normA_p_sq, double *old_err_sq, double *alpha)
{
#ifdef USE_MPI
	double loc_normA_p_sq = *normA_p_sq;
	MPI_Allreduce(&loc_normA_p_sq, normA_p_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

	*alpha = *err_sq / *normA_p_sq ;
	*old_err_sq = *err_sq;

	log_err(SHOW_TASKINFO, "Computing alpha finished : normA_p_sq = %+e ; err_sq = %e ; alpha = %e\n", *normA_p_sq, *err_sq, *alpha);

	// last consumer of these values : let's 0 them so the scalar product doesn't need to
	*err_sq = 0.0;
	*normA_p_sq = 0.0;
}

static inline void swap(double **v, double **w)
{
	double *swap = *v;
	*v = *w;
	*w = swap;
}

void solve_cg(const Matrix *A, const double *b, double *iterate, double convergence_thres)
{
	// do some memory allocations
	double norm_b, thres_sq;
	int r = -1, do_update_gradient = 0;
	double normA_p_sq = 0.0, err_sq = 0.0, old_err_sq = INFINITY, alpha = 0.0, beta = 0.0;

#ifndef USE_MPI
	double *p, *old_p, *Ap, *gradient, *Aiterate;
	p        = (double*)calloc(A->n, sizeof(double));
	old_p    = (double*)calloc(A->n, sizeof(double));
	Ap       = (double*)calloc(A->n, sizeof(double));
	gradient = (double*)calloc(A->n, sizeof(double));
	Aiterate = (double*)calloc(A->n, sizeof(double));
#else
	double *it_glob, *p_glob, *old_p_glob, *p, *old_p, *Ap, *gradient, *Aiterate;

	it_glob    = iterate;
	p_glob     = (double*)calloc(A->m, sizeof(double));
	old_p_glob = (double*)calloc(A->m, sizeof(double));
	Ap         = (double*)calloc(mpi_zonesize[mpi_rank], sizeof(double));
	gradient   = (double*)calloc(mpi_zonesize[mpi_rank], sizeof(double));
	Aiterate   = (double*)calloc(mpi_zonesize[mpi_rank], sizeof(double));

	iterate    = localrank_ptr(it_glob);
	p          = localrank_ptr(p_glob);
	old_p      = localrank_ptr(old_p_glob);
	iterate    = localrank_ptr(it_glob);

	// setting up communications for x and p exchanges
	int first_mpix, last_mpix, count_mpix;
	determine_mpi_neighbours(A, 0, mpi_zonesize[mpi_rank], mpi_size, &first_mpix, &last_mpix);
	count_mpix = last_mpix - first_mpix;

	MPI_Request x_req[2*count_mpix], p1_req[2*count_mpix], p2_req[2*count_mpix], *p_req = p1_req;

	setup_exchange_vect(mpi_rank, first_mpix, last_mpix, 1, it_glob,	x_req);
	setup_exchange_vect(mpi_rank, first_mpix, last_mpix, 2, p_glob,		p1_req);
	setup_exchange_vect(mpi_rank, first_mpix, last_mpix, 3, old_p_glob,	p2_req);
#endif

	// some parameters pre-computed, and show some informations
#ifndef USE_MPI
	norm_b = norm(A->n, b);
#else
	// (borrow thres_sq for local norm_b, get global in norm_b)
	thres_sq = norm(mpi_zonesize[mpi_rank], b);
	MPI_Allreduce(&thres_sq, &norm_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	thres_sq = convergence_thres * convergence_thres * norm_b;
	#if VERBOSE >= SHOW_DBGINFO && (!defined PERFORMANCE || defined EXTRAE_EVENTS)
	log_out("Error shown is ||Ax-b||^2, you should plot ||Ax-b||/||b||. (||b||^2 = %e)\n", norm_b);
	#endif

	#if VERBOSE >= SHOW_TASKINFO
	first_p = p;
	#endif

	start_measure();

	// real work starts here

	for(r=0; r<max_it; r++)
	{
		if( --do_update_gradient > 0 )
		{
			update_gradient(gradient, Ap, &alpha);

			norm_task(gradient, &err_sq);

			compute_beta(&err_sq, &old_err_sq, &beta);

#ifndef USE_MPI
			update_iterate(iterate, old_p, &alpha);
#endif
		}
		else
		{
			// our initial guess is always 0, don't bother updating it
			if( r > 0 )
			{
				update_iterate(iterate, old_p, &alpha);
#ifdef USE_MPI
				#pragma omp task inout({it_glob[glob_block_start(blk):glob_block_end(blk)-1], blk=0:nb_glob_blocks-1}) firstprivate(x_req, count_mpix) label(exchange_x) priority(100) no_copy_deps
				{
					//MPI_Allgatherv(MPI_IN_PLACE, 0/*ignored*/, MPI_DOUBLE, it_glob, mpi_zonesize, mpi_zonestart, MPI_DOUBLE, MPI_COMM_WORLD);

					MPI_Startall(2*count_mpix, x_req);
					MPI_Waitall(2*count_mpix, x_req, MPI_STATUSES_IGNORE);
					log_err(SHOW_TASKINFO, "Exchanging x finished\n");
				}
#endif
			}

#ifndef USE_MPI
			recompute_gradient_mvm(A, iterate, Aiterate);
#else
			recompute_gradient_mvm(A, it_glob, Aiterate);
#endif
			recompute_gradient_update(gradient, Aiterate, b);

			norm_task(gradient, &err_sq);

			compute_beta(&err_sq, &old_err_sq, &beta);
		}

		update_p(p, old_p, gradient, &beta);

#ifdef USE_MPI
		// if possible, execute the update iterate really late, e.g. while we exchange p
		// instead of during the load imbalance while reducing err_sq
		if( do_update_gradient > 0 )
			update_iterate(iterate, old_p, &alpha);

		// need to exchange new p with neighbours
		#pragma omp task inout({p_glob[glob_block_start(blk):glob_block_end(blk)-1], blk=0:nb_glob_blocks-1}) firstprivate(p_req, count_mpix) label(exchange_p) priority(100) no_copy_deps
		{
			// Overlapping opportunity for the runtime work
			//MPI_Allgatherv(MPI_IN_PLACE, 0/*ignored*/, MPI_DOUBLE, p_glob, mpi_zonesize, mpi_zonestart, MPI_DOUBLE, MPI_COMM_WORLD);

			MPI_Startall(2*count_mpix, p_req);
			MPI_Waitall(2*count_mpix, p_req, MPI_STATUSES_IGNORE);

			log_err(SHOW_TASKINFO, "Exchanging p[%p] finished (blocks %d to %d)\n", (void*)(localrank_ptr(p_glob)), world_block(0), world_block(nb_blocks - 1));
		}

		compute_Ap(A, p_glob, Ap);
#else
		compute_Ap(A, p, Ap);
#endif

		scalar_product_task(p, Ap, &normA_p_sq);

#ifndef USE_MPI
		// when reaching this point, all tasks of loop should be created.
		// then waiting starts: should be released halfway through the loop.
		// We want this to be after alpha on normal iterations, after AxIt in recompute iterations
		if( !do_update_gradient )
		{
			#pragma omp taskwait on(old_err_sq) //, {iterate[get_block_start(i):get_block_end(i)-1], i=0:nb_blocks})
		}
		else
		{
			#pragma omp taskwait on(old_err_sq)
		}
#else
		#pragma omp taskwait on(alpha)
#endif
		// swapping p's so we reduce pressure on the execution of the update iterate tasks
		// now output-dependencies are not conflicting with the next iteration but the one after
		{
			swap(&p, &old_p);
#ifdef USE_MPI
			swap(&p_glob, &old_p_glob);
			p_req = (p_req == p1_req) ? p2_req : p1_req;
#endif

			if( r > 0 )
				log_convergence(r-1, old_err_sq);

			if( old_err_sq <= thres_sq )
				break;

			if( do_update_gradient <= 0 )
				do_update_gradient = RECOMPUTE_GRADIENT_FREQ;
		}

		// if we will recompute the gradient, prepare to listen for incoming iterate exchanges in compute_alpha
		compute_alpha(&err_sq, &normA_p_sq, &old_err_sq, &alpha);
	}

	#pragma omp taskwait
	// end of the math, showing infos
	stop_measure();
	log_convergence(r, old_err_sq);


	log_out("CG method finished iterations:%d with error:%e\n", r, sqrt((err_sq==0.0?old_err_sq:err_sq)/norm_b));

#ifdef USE_MPI
	// This is after solving, to be able to compute the verification later on
	//MPI_Allgatherv(MPI_IN_PLACE, 0/*ignored*/, MPI_DOUBLE, it_glob, mpi_zonesize, mpi_zonestart, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Startall(2*count_mpix, x_req);
	MPI_Waitall(2*count_mpix, x_req, MPI_STATUSES_IGNORE);

	for(r=0; r<2*count_mpix; r++)
	{
		MPI_Request_free(x_req+r);
		MPI_Request_free(p1_req+r);
		MPI_Request_free(p2_req+r);
	}

	free(p_glob);
	free(old_p_glob);
#else
	free(p);
	free(old_p);
#endif
	free(Ap);
	free(gradient);
	free(Aiterate);
}

