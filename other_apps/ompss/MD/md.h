#ifndef MD_H
#define MD_H

#ifdef UNUSED
#elif defined(__GNUC__)
	#define UNUSED(x) x __attribute__((unused))
#elif defined(__LCLINT__)
	#define UNUSED(x) /*@unused@*/ x
#else
	#define UNUSED(x) x
#endif

#ifndef DEBUG_PRINT
#define DEBUG_PRINT 0
#endif


typedef struct { double x, y, z; } vect_t;
typedef struct { int x, y, z; } coord_t;

typedef enum {FIXED = 0, PERIODIC = 1} periodicity;

// arrays of atom characteristics
extern vect_t *x, *v, *f;
extern double *m, *ke, *pe;
extern long *id;
// arrays of box characteristics, there are box_count+1 where the last contains lost atoms
long *sb, *nb, box_count;


typedef struct
{
	//space and spatial decomposition
	vect_t lo, hi, size;
	coord_t box_count, period;
	vect_t inv_box_size;
	// containers for atoms and their informations
	long number_atoms, lost_atoms;
} simulation_space;

void allocate(simulation_space *space, double cutoff);

#endif
