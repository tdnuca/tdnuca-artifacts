#include "cg_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

extern int cgver;

extern void *Ahb;
extern void *Acsr;
extern int n;
extern int bm;
extern int cgit;
extern double prec;
extern int correction;
extern int iter_fuse;
extern int rep;
extern char *aname;
extern char *rhsfname;

extern double orth_fac;
extern int cglog_level;
extern int cg_ver;


int cg_config(int argc, char *argv[]) 
{

	if ( argc < 9 ) {
		fprintf(stderr, "use: %s [bm] [it] [precision] [correction] [iter_fuse] [rep] [orth_fac] [HB_MAT] [FULL?] [CG_VER] [cglog] [B]\n", argv[0]);
		return 1;
	}

	rhsfname = NULL;
	cglog_level = 1;

	/* Read parameters */
	bm = atoi(argv[1]);
	cgit = atoi(argv[2]);
	prec = atof(argv[3]);
	correction = atoi(argv[4]);
	iter_fuse = atoi(argv[5]);
	rep = atoi(argv[6]);
	orth_fac = atof(argv[7]);
	aname = argv[8];
	int is_full = atoi(argv[9]);
	cg_ver = atoi(argv[10]);
	if ( argc > 11 ) {
		cglog_level = atoi(argv[11]);
		if ( argc > 12 ) {
			rhsfname = argv[12];
		}
	}

	hbmat_t *tmp = malloc(sizeof(hbmat_t));
	hb_reset(tmp);
	hb_read_double(aname, &(tmp->m), &(tmp->n), &(tmp->elemc), &(tmp->vptr), &(tmp->vpos), (double **)&(tmp->vval));
//	readHB_newmat_double(aname, &(tmp->m), &(tmp->n), &(tmp->elemc), &(tmp->vptr), &(tmp->vpos), (double **)&(tmp->vval));
	tmp->b = bm;
	one2zero(tmp);
	Acsr = tmp;
	if ( !is_full ) {
		Ahb = malloc(sizeof(hbmat_t));
		hb_sym_expand(Ahb, tmp);
		hb_free(tmp);
		Acsr = Ahb;
	}
	n = ((hbmat_t*)Acsr)->m;

	printf("n %d bm %d cgit %d prec %E correction %d iter_fuse %d rep %d orth %E CG_VER %d A %s B %s loglevel %d\n", n, bm, cgit, prec, correction, iter_fuse, rep, orth_fac, cg_ver, aname, rhsfname, cglog_level);

	return 0;
}
