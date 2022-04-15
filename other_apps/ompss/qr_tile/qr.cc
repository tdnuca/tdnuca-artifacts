#include "qr-tile.hh"
#include "misc.hh"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
void usage(char *s)
{
  fprintf(stderr,
	  "Usage: %s <block size> <input file> [output_R output_Q]\nor:\n", s);

  fprintf(stderr,
	  "       %s <block size> <rows> <columns> [output_R output_Q]\n\n", s);
}

int main(int argc, char **argv)
{
  std::ifstream infile;

  if (argc < 3) {
    usage(argv[0]);
    return 1;
  }
  
  // check input file

  Integer BS = atoi(argv[1]);
  if (BS <= 0) {
    fprintf(stderr, 
	    "\nError: Block size must be positive.\n\n");
    usage(argv[0]);
    return 1;
  }


  Integer M, N;

  QR *qr;

  infile.open(argv[2]);
  if (infile.fail()) {
    if (argc < 4) {
      fprintf(stderr, 
	      "\nError: Cannot open \"%s\" for reading, or you have not given the number of columns.\n\n", argv[2]);

      usage(argv[0]);
      return 1;
    }
    M = atol(argv[2]);
    N = atol(argv[3]);

    if (M <= 0 or N <= 0) {
      fprintf(stderr, 
	      "\nError: The number of rows and columns must be positive.\n\n");
      usage(argv[0]);
      return 1;
    }

    qr = new QR(M, N, BS);
    qr->clearMatrix();
    // #pragma omp taskwait
    srand(time(NULL));
    
    Integer NB = qr->A.NB;
    for (Integer I = 0; I < qr->A.NB; ++I) {
      for (Integer J = 0; J < qr->A.MB; ++J) {
	double *m = qr->A(I, J);
#pragma omp task label(init_task) inout([BS*BS]m) firstprivate(I, J, BS, NB)
	{
	  struct drand48_data buf;
	  srand48_r(time(NULL)+I*NB+J, &buf);
	  for (int i = 0; i < BS; ++i) {
	    for (int j = 0; j < BS; ++j) {
	      drand48_r(&buf, &m[i*BS+j]);
	    }
	  }
	}
      }
    }
  }
  else {
    infile >> M >> N;
    qr = new QR(M, N, BS);
    qr->clearMatrix();
#pragma omp taskwait
    for (Integer i = 0; i < M; ++i) {
      for (Integer j = 0; j < N; ++j) {
	infile >> qr->A.elem(i, j);
      }
    }
    infile.close();
  }

  // unsigned np = param.resolution;
  long runtime = wtime();
  qr->run();
  
  runtime = wtime() - runtime;

  fprintf(stdout, "Time: %ld\n", runtime);

  fprintf(stdout, "LAPACK block size: %d\n", qr->LAPACK_bs);

  double flop = 2L*N*N*(M-N/3.0);
  fprintf(stdout, "\t%3.3f GFlop => %6.2f MFlop/s (standard QR)\n", 
	  flop/1000000000.0,
	  flop/runtime);
  fflush(stdout);

  flop = 2L*N*N*(M-N/3.0)*(1L+qr->LAPACK_bs*0.25/BS);
  fprintf(stdout, "\t%3.3f GFlop => %6.2f MFlop/s (true calculations)\n", 
	  flop/1000000000.0,
	  flop/runtime);
  fflush(stdout);
  
  if (argc >= 5) {
    qr->computeQ();
    std::ofstream outfile;
    outfile.open(argv[3]);
    outfile.setf(std::ios::right);
    outfile.precision(15);
    outfile.width(20);
    if (outfile.fail()) {
      fprintf(stderr, 
	      "\nError: Cannot open \"%s\" for writing.\n\n", argv[3]);
      delete qr;
      return 2;
    }
    
    for (Integer i = 0; i < M; ++i) {
      for (Integer j = 0; j < std::min(i, N); ++j) {
	outfile << 0.0 << ' ';
      }
      for (Integer j = i; j < N; ++j) {
	outfile << qr->A.elem(i, j) << ' ';
      }

      outfile << std::endl;
    }

    outfile.close();

#pragma omp taskwait
    outfile.open(argv[4]);
    outfile.setf(std::ios::right);
    outfile.precision(15);
    outfile.width(20);
    if (outfile.fail()) {
      fprintf(stderr, 
	      "\nError: Cannot open \"%s\" for writing.\n\n", argv[4]);

      delete qr;
      return 2;
    }
    
    for (Integer i = 0; i < qr->Q->M; ++i) {
      for (Integer j = 0; j < qr->Q->N; ++j) {
	outfile << qr->Q->elem(i, j) << ' ';
      }

      outfile << std::endl;
    }

    outfile.close();
  }
  delete qr;
}
