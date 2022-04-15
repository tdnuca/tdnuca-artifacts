#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <err.h>

#include <stdint.h>
#include <math.h>
#include <complex.h>

#include "fft.h"
#include "seq.h"


static struct timeval start_time, stop_time;

static inline void start_measure()
{
	gettimeofday(&start_time, NULL);
}


static inline void stop_measure()
{
	gettimeofday(&stop_time, NULL);
	long long time = 1000000 * (stop_time.tv_sec - start_time.tv_sec) + (stop_time.tv_usec - start_time.tv_usec);
	printf("solve_time:%lld\n", time);
}


int open_numpype()
{
	int fd[2];
	if (pipe(fd))
		err(1, "pipe failed");
	int p = fork();

	if (p < 0)
		err(1, "fork failed");

	else if (p == 0)
	{
		// child
		close(fd[1]);
		dup2(fd[0], STDIN_FILENO);
		close(fd[0]);
		char *argv[] = {"python3", NULL};
		execvp(argv[0], argv);
	}

	// parent
	close(fd[0]);
	return fd[1];
}


double lehmer_random(double *seed)
{
	double d2 = .2147483647e10;
	*seed = fmod(16807. * (*seed), d2);
	return (*seed - 1.0) / (d2 - 1.0);
}


void diff_space_to_seq(uint64_t n, double complex *a, double complex *b)
{
	double max_diff = 0, avg_diff = 0;
	for (size_t i = 0; i < n; i++)
	{
		double diff = cabs(b[i] - a[i]);
		if (max_diff < diff)
			max_diff = diff;
		avg_diff += diff / n;
	}

	printf("Over %lu elements, max diff %g avg diff %g\n", n, max_diff, avg_diff);
}


double complex * get_random_matrix(uint64_t n, double seed)
{
	double complex *M = malloc(n * sizeof(*M));

	for (uint64_t i = 0; i < n; i++)
		M[i] = lehmer_random(&seed) + I * lehmer_random(&seed);

	return M;
}


void print_matrix_py(int fd, int ndims, uint64_t dims[ndims], double complex *M)
{
	dprintf(fd, "[");
	if (ndims == 1)
	{
		// NB python tolerates trailing commas just fine
		for (uint64_t j = 0; j < dims[0]; j++)
			dprintf(fd, "%.*g%+.*gj,", __DBL_DECIMAL_DIG__, creal(M[j]), __DBL_DECIMAL_DIG__, cimag(M[j]));
	}
	else if (ndims > 1)
	{
		uint64_t dimsize = 1;
		for (int i = 1; i < ndims; i++)
			dimsize *= dims[i];

		for (uint64_t i = 0; i < dims[0]; i++)
		{
			print_matrix_py(fd, ndims - 1, dims + 1, M + i * dimsize);
			dprintf(fd, ",");
		}
	}
	dprintf(fd, "]");
}


void usage(int exitval, char *argv0)
{
	printf("Usage: %s [-p] [-r] nx [ny [nz]]\n"
	       "Where the size of the d-dimension matrix to transform is nx * ny * nz. Need to be powers of 2.\n"
		   "  -p   : check results using python (numpy.fft.ffn())\n"
		   "  -r   : check results using seq\n",
	       argv0);
	exit(exitval);
}

int main(int argc, char *argv[])
{
	int numpy_check = 0, seq_check = 0, ndims = 0, opt;
	double seed = 564321;

	while ((opt = getopt(argc, argv, "rps:")) != -1)
	{
		if (opt == 'p')
			numpy_check = 1;
		else if (opt == 'r')
			seq_check = 1;
		else if (opt == 's')
			seed = strtod(optarg, NULL);
	}

	ndims = argc - optind;
	if (ndims < 1 || ndims > 3)
		usage(0, argv[0]);


	uint64_t dims[ndims], total = 1;
	for (int i = 0; i < ndims; i++)
	{
		dims[i] = strtoull(argv[i + optind], NULL, 0);
		if (!dims[i] || dims[i] != (1ULL << (ffs(dims[i]) - 1)))
			usage(1, argv[0]);
		total *= dims[i];
	}


	complex double *space = get_random_matrix(total, seed), *seq_space = NULL;

	int fd;
	if (numpy_check)
	{
		fd = open_numpype();
		dprintf(fd, "import numpy as np\n");
		dprintf(fd, "a = np.array(");
		print_matrix_py(fd, ndims, dims, space);
		dprintf(fd, ")\n");
	}

	if (seq_check)
		seq_space = memcpy(malloc(total * sizeof(*seq_space)), space, total * sizeof(*space));


	start_measure();

	if (ndims == 3)
		fft3D(dims[0], dims[1], dims[2], (complex double (*)[dims[1]][dims[2]])space, 1);
	if (ndims == 2)
		fft2D(dims[0], dims[1], (complex double (*)[dims[1]])space, 1);
	if (ndims == 1)
		fft(dims[0], space, 1);

	stop_measure();


	if (seq_check)
	{
		// compare against timing & results of sequential implementation
		start_measure();

		if (ndims == 3)
			seq_fft3D(dims[0], dims[1], dims[2], (complex double (*)[dims[1]][dims[2]])seq_space, 1);
		if (ndims == 2)
			seq_fft2D(dims[0], dims[1], (complex double (*)[dims[1]])seq_space, 1);
		if (ndims == 1)
			seq_fft(dims[0], seq_space, 1);

		stop_measure();
		diff_space_to_seq(total, space, seq_space);
	}

	if (numpy_check)
	{
		dprintf(fd, "b = np.array(");
		print_matrix_py(fd, ndims, dims, space);
		dprintf(fd, ")\n");
		dprintf(fd, "print(\"max error is\", np.absolute(np.fft.fftn(a) - b).max())\n");

		close(fd);
		wait(NULL);
	}

	free(space);
	free(seq_space);

	return 0;
}
