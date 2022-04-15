#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <float.h>

#include "md.h"
#include "io_lammps.h"

vect_t *x, *v, *f;
double *m, *ke, *pe;
long *id, *sb, *nb, box_count;

void sum_vect(vect_t *out, vect_t *in)
{
	out->x += in->x;
	out->y += in->y;
	out->z += in->z;
}

void max_vect(vect_t *out, vect_t *in)
{
	double o_sq = out->x * out->x + out->y * out->y + out->z * out->z;
	double i_sq = in->x * in->x + in->y * in->y + in->z * in->z;
	if (o_sq < i_sq)
		*out = *in;
}

#pragma omp declare reduction(+:   vect_t: sum_vect(&omp_out, &omp_in)) initializer(omp_priv = {0})
#pragma omp declare reduction(max: vect_t: max_vect(&omp_out, &omp_in)) initializer(omp_priv = {0})

static inline int sgn(int x) { return (x > 0) - (x < 0); }

void allocate(simulation_space *space, double cutoff)
{
	m = calloc(space->number_atoms, sizeof(*m));
	x = calloc(space->number_atoms, sizeof(*x));
	v = calloc(space->number_atoms, sizeof(*v));
	f = calloc(space->number_atoms, sizeof(*f));
	id = calloc(space->number_atoms, sizeof(*id));
	ke = calloc(space->number_atoms, sizeof(*ke));
	pe = calloc(space->number_atoms, sizeof(*pe));

	space->size = (vect_t)
	{
		.x = space->hi.x - space->lo.x,
		.y = space->hi.y - space->lo.y,
		.z = space->hi.z - space->lo.z
	};

	space->box_count = (coord_t)
	{
		.x = floor(space->size.x / cutoff),
		.y = floor(space->size.y / cutoff),
		.z = floor(space->size.z / cutoff)
	};

	box_count = space->box_count.x * space->box_count.y * space->box_count.z;
	sb = calloc(box_count + 1, sizeof(*sb));
	nb = calloc(box_count + 1, sizeof(*nb));

	space->inv_box_size = (vect_t)
	{
		.x = (double)space->box_count.x / space->size.x,
		.y = (double)space->box_count.y / space->size.y,
		.z = (double)space->box_count.z / space->size.z
	};

	// safety measures. If box size 10, box count 3, an atom at -eps might be moved to 10-eps -> 10
	// then 10*(1/3) needs to be < 3
	int pow;
	while (space->size.x * space->inv_box_size.x >= space->box_count.x)
	{
		frexp(space->inv_box_size.x, &pow);
		space->inv_box_size.x -= ldexp(DBL_EPSILON, pow);
	}
	while (space->size.y * space->inv_box_size.y >= space->box_count.y)
	{
		frexp(space->inv_box_size.y, &pow);
		space->inv_box_size.y -= ldexp(DBL_EPSILON, pow);
	}
	while (space->size.z * space->inv_box_size.z >= space->box_count.z)
	{
		frexp(space->inv_box_size.z, &pow);
		space->inv_box_size.z -= ldexp(DBL_EPSILON, pow);
	}

	assert(space->size.x * space->inv_box_size.x < space->box_count.x);
	assert(space->size.y * space->inv_box_size.y < space->box_count.y);
	assert(space->size.z * space->inv_box_size.z < space->box_count.z);
	if (1 / space->inv_box_size.x < cutoff || 1 / space->inv_box_size.y < cutoff || 1 / space->inv_box_size.z < cutoff)
		fprintf(stderr, "WARNING your box size is [%g,%g,%g] but cutoff is %g\n", 1 / space->inv_box_size.x, 1 / space->inv_box_size.y, 1 / space->inv_box_size.z, cutoff);

}

void deallocate(simulation_space *space)
{
	(void)space;
	free(m);
	free(x);
	free(v);
	free(f);
	free(id);
	free(ke);
	free(pe);
	free(sb);
	free(nb);
}

static inline void single_force(int atom_a, int atom_b, vect_t skew, double cutoff2, double epsilon, double sigma2)
{
	double dx = x[atom_a].x - (x[atom_b].x + skew.x);
	double dy = x[atom_a].y - (x[atom_b].y + skew.y);
	double dz = x[atom_a].z - (x[atom_b].z + skew.z);

	double d2 = dx * dx + dy * dy + dz * dz;
	if (d2 > cutoff2)
		return;

	assert(atom_a != atom_b);

	double r2 = 1.0 / d2;
	double rs2 = r2 * sigma2;
	double rs6 = rs2 * rs2 * rs2 ;

	double norm_fab = 24 * epsilon * rs6 * (1 - 2 * rs6) * r2;
	double potential = 4 * epsilon * rs6 * (rs6 - 1) / 2;

	// force of particle a on particle b
	double fab_x = dx * norm_fab;
	double fab_y = dy * norm_fab;
	double fab_z = dz * norm_fab;

	// NB all of these *would need* to be atomic, if force tasks were called from within a concurrent() task
	f[atom_b].x += fab_x;
	f[atom_b].y += fab_y;
	f[atom_b].z += fab_z;

	f[atom_a].x -= fab_x;
	f[atom_a].y -= fab_y;
	f[atom_a].z -= fab_z;

	pe[atom_b] += potential;
	pe[atom_a] += potential;

	if (DEBUG_PRINT)
		fprintf(stderr, "d2 = %g r2 = %g sigma2 = %g rs2 = %g rs6 = %g pij = %g\n", d2, r2, sigma2, rs2, rs6, potential * 2);
	if (DEBUG_PRINT)
		fprintf(stderr,"Norm of the force between %ld (%g,%g,%g) and %ld (%g,%g,%g) is %g : force is (%g,%g,%g)\nNow f(%ld) is (%g,%g,%g) and f(%ld) is (%g,%g,%g)\n",
				id[atom_a], x[atom_a].x, x[atom_a].y, x[atom_a].z, id[atom_b], x[atom_b].x + skew.x, x[atom_b].y + skew.y, x[atom_b].z + skew.z,
				norm_fab * sqrt(d2), fab_x, fab_y, fab_z, id[atom_a], f[atom_a].x, f[atom_a].y, f[atom_a].z, id[atom_b], f[atom_b].x, f[atom_b].y, f[atom_b].z);
}


void self_force(long box, double cutoff2,  double epsilon, double sigma2)
{
	memset(f + sb[box], 0, nb[box] * sizeof(*f));
	memset(pe + sb[box], 0, nb[box] * sizeof(*pe));

	vect_t noskew = {0};
	for (int i = 0; i < nb[box]; i++)
		for (int j = i + 1; j < nb[box]; j++)
			single_force(sb[box] + i, sb[box] + j, noskew, cutoff2, epsilon, sigma2);
}

void force(long box_a, long box_b, vect_t skew, double cutoff2,  double epsilon, double sigma2)
{
	for (int i = 0; i < nb[box_a]; i++)
		for (int j = 0; j < nb[box_b]; j++)
			single_force(sb[box_a] + i, sb[box_b] + j, skew, cutoff2, epsilon, sigma2);
}

static inline
long box_id(vect_t *pos, simulation_space *space)
{
	// wrap_d is 0 if d is a non-periodic dimension (no wrapping) or coordinate inside bounds, 1 if too low, -1 if too high
	// Arguably the speed should never be such that we escape the simulation space by more than the size of the space
	int wrap_x = space->period.x * ((pos->x < space->lo.x) - (pos->x > space->hi.x));
	int wrap_y = space->period.y * ((pos->y < space->lo.y) - (pos->y > space->hi.y));
	int wrap_z = space->period.z * ((pos->z < space->lo.z) - (pos->z > space->hi.z));

	pos->x += wrap_x * space->size.x;
	pos->y += wrap_y * space->size.y;
	pos->z += wrap_z * space->size.z;

	long x = floor((pos->x - space->lo.x) * space->inv_box_size.x);
	long y = floor((pos->y - space->lo.y) * space->inv_box_size.y);
	long z = floor((pos->z - space->lo.z) * space->inv_box_size.z);

	/* max value +1 to indicate out of bounds */
	if (x < 0 || y < 0 || z < 0 || x >= space->box_count.x || y >= space->box_count.y || z >= space->box_count.z)
		return space->box_count.x * space->box_count.y * space->box_count.z;

	return (x * space->box_count.y + y) * space->box_count.z + z;
}

void integrate(long box, const double dtime)
{
	const double v_dtime = dtime / 2;

	for (int i = sb[box]; i < sb[box] + nb[box]; i++)
	{
		v[i].x += v_dtime * f[i].x / m[i];
		v[i].y += v_dtime * f[i].y / m[i];
		v[i].z += v_dtime * f[i].z / m[i];

		x[i].x += dtime * v[i].x;
		x[i].y += dtime * v[i].y;
		x[i].z += dtime * v[i].z;

		if (DEBUG_PRINT)
			fprintf(stderr, "Integrated position of i %ld to\tx[%g,%g,%g]\tv[%g,%g,%g]\tf[%g,%g,%g]\n",
					id[i], x[i].x, x[i].y, x[i].z, v[i].x, v[i].y, v[i].z, f[i].x, f[i].y, f[i].z);
	}
}

void progress_speed_intstep(long box, double dtime)
{
	const double v_dtime = dtime / 2;

	for (int i = sb[box]; i < sb[box] + nb[box]; i++)
	{
		v[i].x += v_dtime * f[i].x / m[i];
		v[i].y += v_dtime * f[i].y / m[i];
		v[i].z += v_dtime * f[i].z / m[i];

		// kinetic energy, mv^2
		ke[i] = (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z) * m[i] / 2;
	}
}

static inline
void swap_atoms(int i, int j)
{
	// save atom at i, to fill position with atom that goes here
	vect_t t_x = x[i];
	vect_t t_v = v[i];
	double t_m = m[i];
	int t_id = id[i];

	x[i] = x[j];
	v[i] = v[j];
	m[i] = m[j];
	id[i] = id[j];

	x[j] = t_x;
	v[j] = t_v;
	m[j] = t_m;
	id[j] = t_id;
}

void sort_atoms(simulation_space *space, int init)
{
	int *order = malloc(space->number_atoms * sizeof(*order));
	memset(nb, 0, (box_count + 1) * sizeof(*nb));

	// A ϴ(N) bucket-ish sort: assign boxes to atoms & get box sizes
    for(int i = 0; i < space->number_atoms; i++)
	{
		int box = box_id(x + i, space);
		nb[box]++;
		order[i] = box;

		if (DEBUG_PRINT && init)
			fprintf(stderr, "Atom %ld [%.2g,%.2g,%.2g] added to box [%d, %d, %d]\n", id[i], x[i].x, x[i].y, x[i].z,
					box / (space->box_count.z * space->box_count.y), (box / space->box_count.z) % space->box_count.y, box % space->box_count.z);
	}

	// in order of reading for sorting stability, assign each atom to a position in its box
	long boxpos[box_count + 1];
	boxpos[0] = sb[0] = 0;
	for (int i = 1; i < box_count + 1; i++)
		boxpos[i] = sb[i] = sb[i - 1] + nb[i - 1];

    for(int i = 0; i < space->number_atoms; i++)
		order[i] = boxpos[order[i]]++;

	// once the positions are assigned, move atoms based on ordering in order[]
	for (long i = 0; i < space->number_atoms; )
	{
		if (order[i] == i)
			i++;
		else
		{
			int dest = order[i];
			swap_atoms(i, dest);

			order[i] = order[dest];
			order[dest] = dest;
		}
	}

	free(order);
}

int get_half_neighbours(int *neigh, coord_t *wrap, coord_t pos, coord_t box_count, coord_t period)
{
	int n_neigh = 0, i, j, k, ii, jj, kk;
	coord_t wrapping;

	// not -1, removing 9 half-neighbours
	for (i = pos.x; i < pos.x + 2; i++)
	{
		if (i == -1 || i == box_count.x)
		{
			if (!period.x)
				continue;

			wrapping.x = sgn(i);
			ii = (i + box_count.x) % box_count.x;
		}
		else
		{
			wrapping.x = 0;
			ii = i;
		}

		for (j = pos.y - 1; j < pos.y + 2; j++)
		{
			// not -1 on i==x, removing 3 half neighbours
			if (i == pos.x && j + 1 == pos.y)
				continue;
			else if (j == -1 || j == box_count.y)
			{
				if (!period.y)
					continue;

				wrapping.y = sgn(j);
				jj = (j + box_count.y) % box_count.y;
			}
			else
			{
				wrapping.y = 0;
				jj = j;
			}

			for (k = pos.z - 1; k < pos.z + 2; k++)
			{
				if (k == -1 || k == box_count.z)
				{
					if (!period.z)
						continue;

					wrapping.z = sgn(k);
					kk = (k + box_count.z) % box_count.z;
				}
				else
				{
					wrapping.z = 0;
					kk = k;
				}

				// yay ! But only for points not central, or the last half-neighbour
				if (i == pos.x && j == pos.y && (k == pos.z || k + 1 == pos.z))
					continue;

				wrap[n_neigh] = wrapping;
				neigh[n_neigh] = (ii * box_count.y + jj) * box_count.z + kk;
				n_neigh++;
			}
		}
	}

	if (DEBUG_PRINT)
	{
		fprintf(stderr, "box [%d,%d,%d] has half-neighbours ", pos.x, pos.y, pos.z);
		for (i = 0; i < n_neigh; i++)
			fprintf(stderr, "[%d,%d,%d] ", neigh[i] / (box_count.z * box_count.y), (neigh[i] / box_count.z) % box_count.y, neigh[i] % box_count.z);
		fprintf(stderr, "\n");
	}

	return n_neigh;
}

void usage(const char *argv0)
{
	fprintf(stderr, "Usage : %s [-dt=0.005] [-epsilon=1.0] [-sigma=1.0] [-cutoff=2.5] [-it=100] [-write_freq=100] [-write_file=out.dump] /path/to/file.dump\n", argv0);
	exit(-1);
}

int main(int argc, char *argv[])
{
	simulation_space space;
	// few things to read from command line, put in default values

	space.lost_atoms = 0;
	double sigma = 1, epsilon = 1;
	double dtime = 0.001, cutoff = 2.5; // LJ
	int time = 0, iterations = 100, write_freq = iterations;
	char *filename = NULL;
	FILE *out = stdout;

	for (int n = 1; n < argc; n++)
	{
		if (strstr(argv[n], "-dt=") == argv[n])
			dtime = strtod(strchr(argv[n], '=') + 1, NULL);
		else if (strstr(argv[n], "-epsilon=") == argv[n])
			epsilon = strtod(strchr(argv[n], '=') + 1, NULL);
		else if (strstr(argv[n], "-sigma=") == argv[n])
			sigma = strtod(strchr(argv[n], '=') + 1, NULL);
		else if (strstr(argv[n], "-cutoff=") == argv[n])
			cutoff = strtod(strchr(argv[n], '=') + 1, NULL);
		else if (strstr(argv[n], "-it=") == argv[n])
			iterations = atoi(strchr(argv[n], '=') + 1);
		else if (strstr(argv[n], "-write_freq=") == argv[n])
			write_freq = atoi(strchr(argv[n], '=') + 1);
		else if (strstr(argv[n], "-write_file") == argv[n])
			out = fopen(strchr(argv[n], '=') + 1, "w");
		else if (strchr(argv[n], '=') == NULL)
			filename = argv[n];
		else
		{
			fprintf(stderr, "Unrecognized option %s\n", argv[n]);
			usage(argv[0]);
		}
	}

	if (filename == NULL)
	{
		fprintf(stderr, "I don't know what the atoms I should simulate are. Feed me a LAMMPS dump file !\n");
		usage(argv[0]);
	}

	FILE *h = read_sim_metadata(filename, &space, &time);

	// allocate & load atoms
	allocate(&space, cutoff);
	read_atoms(h, space.number_atoms);

	// some more initializations
	double cutoff2 = cutoff * cutoff, sigma2 = sigma * sigma, kinetic_energy, potential_energy;

	// sum up parameters, and off to real work
	fprintf(stderr, "Using dt = %g, cutoff = %g, nb. boxes in space [%d, %d, %d] iterations %d, dumping every %d\n", dtime, cutoff, space.box_count.x, space.box_count.y, space.box_count.z, iterations, write_freq);

	// neighbour boxes
	int half_neigh[box_count][13], n_half[box_count];
	coord_t wrap[box_count][26];

	for (int i = 0, n = 0; i < space.box_count.x; i++)
		for (int j = 0; j < space.box_count.y; j++)
			for (int k = 0; k < space.box_count.z; k++, n++)
				n_half[n] = get_half_neighbours(half_neigh[n], wrap[n], (coord_t) {.x = i, .y = j, .z = k}, space.box_count, space.period);


	sort_atoms(&space, 1);

	int t;
	for (t = 0; t <= iterations; t++)
	{
		fprintf(stderr, "ITERATION %d\n", t);

		for (int n = 0; n < box_count; n++)
			if (nb[n])
				// forces of particles in the box on each other
				#pragma omp task out(f[sb[n];nb[n]], pe[sb[n];nb[n]]) in(x[sb[n];nb[n]]) firstprivate(n) label(self_force)
				self_force(n, cutoff2, epsilon, sigma2);

		for (int n = 0; n < box_count; n++)
			for (int l = 0; l < n_half[n]; l++)
			{
				vect_t skew = (vect_t) {.x = wrap[n][l].x * space.size.x, .y = wrap[n][l].y * space.size.y, .z = wrap[n][l].z * space.size.z};
				int h = half_neigh[n][l];

				if (!nb[n] || !nb[h])
					continue;

				// Nanox reductions use OpenMP syntax for some reason, so spoon-feed them a [size]pointer format
				vect_t *UNUSED(fn)  = f  + sb[n], *UNUSED(fm) = f + sb[h];
				double *UNUSED(pen) = pe + sb[n], *UNUSED(pem) = pe + sb[h];

				#pragma omp task reduction(+: [nb[n]]fn, [nb[h]]fm, [nb[n]]pen, [nb[h]]pem) in(x[sb[n];nb[n]], x[sb[h];nb[h]]) label(force)
				force(n, h, skew, cutoff2, epsilon, sigma2);
			}


		for (int n = 0; n < box_count; n++)
			if (nb[n])
				#pragma omp task inout(v[sb[n];nb[n]]) out(ke[sb[n];nb[n]]) in(f[sb[n];nb[n]], m[sb[n];nb[n]]) label(integrate_speed)
				progress_speed_intstep(n, t ? dtime : 0.0);

		if (t == iterations || t % write_freq == 0)
		{
			kinetic_energy = 0.;
			potential_energy = 0.;

			for (int n = 0; n < box_count; n++)
				#pragma omp task reduction(+: kinetic_energy, potential_energy) in(ke[sb[n];nb[n]], pe[sb[n];nb[n]]) label(sum_energy)
				for (int i = sb[n]; i < sb[n] + nb[n]; i++)
				{
					kinetic_energy += ke[i];
					potential_energy += pe[i];
				}

			#pragma omp taskwait
			write_out(out, &space, time + t, kinetic_energy, potential_energy);
		}

		for (int n = 0; n < box_count; n++)
			if (nb[n])
				#pragma omp task inout(x[sb[n];nb[n]], v[sb[n];nb[n]]) in(f[sb[n];nb[n]], m[sb[n];nb[n]]) label(integrate_position)
				integrate(n, dtime);

		#pragma omp taskwait
		sort_atoms(&space, 0);
	}

	// just for the sake of separating iterations
	#pragma omp taskwait

	// if not periodic or some way too high velocity atoms escaped, warn the user
	if ((!space.period.x || !space.period.y || !space.period.z) || space.lost_atoms)
		fprintf(stderr, "You have lost %ld atoms during this simulation.\n", space.lost_atoms);

	deallocate(&space);
	if (out != stdout)
		fclose(out);

	return 0;
}

