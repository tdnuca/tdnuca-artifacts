#ifndef READ_LAMMPS_H
#define READ_LAMMPS_H

#include "md.h"

FILE *read_sim_metadata(char *filename, simulation_space *space, int *time);
void read_atoms(FILE *handle, int n_atoms);

void write_out(FILE *handle, simulation_space *space, int time, double kinetic, double potential);
void write_file(char *filename, simulation_space *space, int time, double kinetic, double potential);

#endif

