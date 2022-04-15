#ifndef __CG_SETUP_H__
#define __CG_SETUP_H__


#include "cg_aux.h"
#include "vector.h"
#include "csparse.h"


int cg_setup(int n, int bm, double **x, double **rhs);


#endif // __CG_SETUP_H__
