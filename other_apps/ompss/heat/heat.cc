#include "heat.hh"
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#ifdef ENABLE_PARSEC_HOOKS
#include "hooks.h"
#endif

void usage(char *s)
{
  fprintf(stderr,
	  "Usage: %s <number of blocks> <input file> [result file]\n\n", s);
}

int main(int argc, char **argv)
{
  unsigned long iter;
  std::ifstream infile;
  std::FILE *resfile;
  std::string resfilename;

  algoparam_t param;
  
  double runtime, flop;
  double residual = 0.0;

  double runtime2 = wtime();

  if (argc < 3) {
    usage(argv[0]);
    return 1;
  }

  

  // check input file
  infile.open(argv[2]);
  if (infile.fail()) {
    fprintf(stderr, 
	    "\nError: Cannot open \"%s\" for reading.\n\n", argv[2]);
    
    usage(argv[0]);
    return 1;
  }

  if (not param.read_input(infile)) {
    infile.close();
    fprintf(stderr, "\nError: Error parsing input file.\n\n");
    usage(argv[0]);
    return 1;
  }
  infile.close();

  if (param.algorithm < 0 or param.algorithm >= 3) {
    std::cerr << "Unknown algorithm selected (" << param.algorithm << ')' << std::endl;
    return 1;
  }

  // param.print_params;

  if (not param.init(atol(argv[1]))) {
    fprintf(stderr, "Error in Solver initialization.\n\n");
    usage(argv[0]);
    return 1;
  }
  
  param.print_params();

  // unsigned np = param.resolution;
  runtime = wtime();

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_begin();
#endif
  
  iter = 1;
  while (true) {
    bool check = (iter%1000)==0 or (iter==param.maxiter);

    switch (param.algorithm) {
    case 0:
      relax_jacobi(*param.m, check, &residual);
      break;
    case 1:
      relax_gauss(*param.m, check, &residual);
      break;
    case 2:
      relax_redblack(*param.m, check, &residual);
      break;
    default:
      std::cerr << "Unknown algorithm selected (" << param.algorithm << ')' << std::endl;
      return 1;
    }
    
    if (check) {
#pragma omp taskwait
      if (residual < 0.00005)
	break;
    }

    ++iter;
    
    if (param.maxiter > 0 and iter > param.maxiter) {
      --iter;
      break;
    }
  }
  // not needed because of `check' usage
  // #pragma omp taskwait

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif

  flop = iter*11.0*param.resolution*param.resolution;
  runtime = wtime() - runtime;
  runtime2 = wtime() - runtime2;

  fprintf(stdout, "Time: %04.3f ", runtime);
  fprintf(stdout, " /-> Time total: %04.3f ", runtime2);
  
  fprintf(stdout, "(%3.3f GFlop => %6.2f MFlop/s)\n", 
	  flop/1000000000.0,
	  flop/runtime/1000000);
  fprintf(stdout, "Convergence to residual=%f: %lu iterations\n", residual, iter);
  
  // for plot...
  if (argc >= 4) {
    resfilename = std::string(argv[3]);
    resfile = std::fopen(resfilename.c_str(), "w");
    if (resfile == 0) {
      fprintf(stderr, 
	      "\nError: Cannot open \"%s\" for writing.\n\n", 
	      resfilename.c_str());
      usage(argv[0]);
      return 1;
    }
    
    param.coarsen();
    param.write_image(resfile);
    std::fclose(resfile);
  }

  return 0;
}
