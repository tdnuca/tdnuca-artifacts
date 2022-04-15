#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>

#include "md.h"

typedef enum lammps_attribute
{
	LAMMPS_FIELD_UNKNOWN = 0, LAMMPS_FIELD_ATOMID = 1, LAMMPS_FIELD_MOLID, LAMMPS_FIELD_TYPE,
	LAMMPS_FIELD_POSX = 3, LAMMPS_FIELD_POSY = 4, LAMMPS_FIELD_POSZ = 5, LAMMPS_FIELD_POSXS,
	LAMMPS_FIELD_POSYS, LAMMPS_FIELD_POSZS, LAMMPS_FIELD_POSXU, LAMMPS_FIELD_POSYU,
	LAMMPS_FIELD_POSZU, LAMMPS_FIELD_POSXSU, LAMMPS_FIELD_POSYSU, LAMMPS_FIELD_POSZSU,
	LAMMPS_FIELD_IMGX, LAMMPS_FIELD_IMGY, LAMMPS_FIELD_IMGZ, LAMMPS_FIELD_VELX = 6,
	LAMMPS_FIELD_VELY = 7, LAMMPS_FIELD_VELZ = 8, LAMMPS_FIELD_FORX, LAMMPS_FIELD_FORY,
	LAMMPS_FIELD_FORZ, LAMMPS_FIELD_CHARGE, LAMMPS_FIELD_RADIUS, LAMMPS_FIELD_DIAMETER,
	LAMMPS_FIELD_ELEMENT, LAMMPS_FIELD_MASS = 2, LAMMPS_FIELD_QUATW, LAMMPS_FIELD_QUATI,
	LAMMPS_FIELD_QUATJ, LAMMPS_FIELD_QUATK, LAMMPS_FIELD_MUX, LAMMPS_FIELD_MUY,
	LAMMPS_FIELD_MUZ, LAMMPS_FIELD_USER0, LAMMPS_FIELD_USER1, LAMMPS_FIELD_USER2,
	LAMMPS_FIELD_USER3, LAMMPS_FIELD_USER4, LAMMPS_FIELD_USER5, LAMMPS_FIELD_USER6,
	LAMMPS_FIELD_USER7, LAMMPS_FIELD_USER8, LAMMPS_FILED_USER9
} l_attr_t;


FILE *read_sim_metadata(char *filename, simulation_space *space, int *time)
{
	// some of this is too verbose but copied from lammps source code
	FILE *handle = fopen(filename, "r");

	if (handle == NULL)
		err(1, "Failed opening file %s", filename);

	char line[250];
	fpos_t start_of_line;

	while (!feof(handle) && fgetpos(handle, &start_of_line) == 0 && fgets(line, sizeof(line), handle))
	{
		if (strstr(line, "ITEM: TIMESTEP"))
		{
			*time = strtol(fgets(line, sizeof(line), handle), NULL, 10);
			if (*time < 0)
				err(3, "Failed to read time from file %s", filename);
		}

		else if (strstr(line, "ITEM: NUMBER OF ATOMS"))
		{
			space->number_atoms = strtol(fgets(line, sizeof(line), handle), NULL, 10);
			if (space->number_atoms <= 0)
				err(4, "Failed to read number of atoms from file %s", filename);
		}

		else if (strstr(line, "ITEM: BOX BOUNDS"))
		{
			char *k = strtok(line + 17, " \t\n\r"); // skip "ITEM: BOX BOUNDS "
			periodicity per[3] = {FIXED};

			for (int i = 0; i < 3 && k != NULL; i++, k = strtok(NULL, " \t\n\r"))
			{
				if (0 == strcmp(k, "ff"))
					per[i] = FIXED;
				else if (0 == strcmp(k, "pp"))
					per[i] = PERIODIC;
				else
					errx(5, "Unexpected periodicity type in file %s: %s", filename, k);
			}

			space->period = (coord_t){.x = per[0], .y = per[1], .z = per[2]};

			char *next;
			space->lo.x = strtod(fgets(line, sizeof(line), handle), &next);
			space->hi.x = strtod(next, NULL);
			space->lo.y = strtod(fgets(line, sizeof(line), handle), &next);
			space->hi.y = strtod(next, NULL);
			space->lo.z = strtod(fgets(line, sizeof(line), handle), &next);
			space->hi.z = strtod(next, NULL);
			space->size = (vect_t) {.x = space->hi.x - space->lo.x, .y = space->hi.y - space->lo.y, .z = space->hi.z - space->lo.z};
		}

		// reached the atoms: rewind this line and stop reading here
		else if (strstr(line, "ITEM: ATOMS"))
		{
			fsetpos(handle, &start_of_line);
			return handle;
		}
	}

	fclose(handle);

	errx(6, "Reached end of file before finding atoms!\n");
	return NULL;
}


void read_atoms(FILE *handle, long n_atoms)
{
	char line[250];

	fgets(line, sizeof(line), handle);
	if (!strstr(line, "ITEM: ATOMS"))
		errx(7, "Could not find start of atoms in file\n");

	// headers for fields that describe atoms
	char *k = strtok(line + 11, " \t\n\r"); // skip "ITEM: ATOMS "
	l_attr_t field_types[50] = {LAMMPS_FIELD_UNKNOWN};
	int number_fields = 0;

	for (int i = 0; k != NULL; i++, k = strtok(NULL, " \t\n\r"))
	{
		if (0 == strcmp(k, "id"))              field_types[i] = LAMMPS_FIELD_ATOMID;
		else if (0 == strcmp(k, "x"))          field_types[i] = LAMMPS_FIELD_POSX;
		else if (0 == strcmp(k, "y"))          field_types[i] = LAMMPS_FIELD_POSY;
		else if (0 == strcmp(k, "z"))          field_types[i] = LAMMPS_FIELD_POSZ;
		else if (0 == strcmp(k, "vx"))         field_types[i] = LAMMPS_FIELD_VELX;
		else if (0 == strcmp(k, "vy"))         field_types[i] = LAMMPS_FIELD_VELY;
		else if (0 == strcmp(k, "vz"))         field_types[i] = LAMMPS_FIELD_VELZ;
		else if (0 == strcmp(k, "mass"))       field_types[i] = LAMMPS_FIELD_MASS;
		//else if (0 == strcmp(k, "mol"))      field_types[i] = LAMMPS_FIELD_MOLID;
		//else if (0 == strcmp(k, "type"))     field_types[i] = LAMMPS_FIELD_TYPE;
		//else if (0 == strcmp(k, "xs"))       field_types[i] = LAMMPS_FIELD_POSXS;
		//else if (0 == strcmp(k, "ys"))       field_types[i] = LAMMPS_FIELD_POSYS;
		//else if (0 == strcmp(k, "zs"))       field_types[i] = LAMMPS_FIELD_POSZS;
		//else if (0 == strcmp(k, "xu"))       field_types[i] = LAMMPS_FIELD_POSXU;
		//else if (0 == strcmp(k, "yu"))       field_types[i] = LAMMPS_FIELD_POSYU;
		//else if (0 == strcmp(k, "zu"))       field_types[i] = LAMMPS_FIELD_POSZU;
		//else if (0 == strcmp(k, "xus"))      field_types[i] = LAMMPS_FIELD_POSXU;
		//else if (0 == strcmp(k, "yus"))      field_types[i] = LAMMPS_FIELD_POSYU;
		//else if (0 == strcmp(k, "zus"))      field_types[i] = LAMMPS_FIELD_POSZU;
		//else if (0 == strcmp(k, "ix"))       field_types[i] = LAMMPS_FIELD_IMGX;
		//else if (0 == strcmp(k, "iy"))       field_types[i] = LAMMPS_FIELD_IMGY;
		//else if (0 == strcmp(k, "iz"))       field_types[i] = LAMMPS_FIELD_IMGZ;
		//else if (0 == strcmp(k, "fx"))       field_types[i] = LAMMPS_FIELD_FORX;
		//else if (0 == strcmp(k, "fy"))       field_types[i] = LAMMPS_FIELD_FORY;
		//else if (0 == strcmp(k, "fz"))       field_types[i] = LAMMPS_FIELD_FORZ;
		//else if (0 == strcmp(k, "q"))        {field_types[i] = LAMMPS_FIELD_CHARGE; *optflags |= MOLFILE_CHARGE;}
		//else if (0 == strcmp(k, "radius"))   {field_types[i] = LAMMPS_FIELD_RADIUS; *optflags |= MOLFILE_RADIUS;}
		//else if (0 == strcmp(k, "diameter")) field_types[i] = LAMMPS_FIELD_RADIUS;
		//else if (0 == strcmp(k, "element"))  field_types[i] = LAMMPS_FIELD_ELEMENT;
		//else if (0 == strcmp(k, "mux"))      field_types[i] = LAMMPS_FIELD_MUX;
		//else if (0 == strcmp(k, "muy"))      field_types[i] = LAMMPS_FIELD_MUY;
		//else if (0 == strcmp(k, "muz"))      field_types[i] = LAMMPS_FIELD_MUZ;

		if (field_types[i] != LAMMPS_FIELD_UNKNOWN)
			number_fields++;
	}

	// read the atoms (one per line)
	long i;
	for (i = 0; !feof(handle) && fgets(line, sizeof(line), handle); i++)
	{
		int used_fields = 0;
		l_attr_t *f = field_types;
		for (char *k = strtok(line, " \t"); k != NULL; k = strtok(NULL, " \t"), f++)
		{
			if (*f == LAMMPS_FIELD_ATOMID) id[i] = strtol(k, NULL, 10);
			else if (*f == LAMMPS_FIELD_POSX) x[i].x = strtod(k, NULL);
			else if (*f == LAMMPS_FIELD_POSY) x[i].y = strtod(k, NULL);
			else if (*f == LAMMPS_FIELD_POSZ) x[i].z = strtod(k, NULL);
			else if (*f == LAMMPS_FIELD_VELX) v[i].x = strtod(k, NULL);
			else if (*f == LAMMPS_FIELD_VELY) v[i].y = strtod(k, NULL);
			else if (*f == LAMMPS_FIELD_VELZ) v[i].z = strtod(k, NULL);
			else if (*f == LAMMPS_FIELD_MASS) m[i] = strtod(k, NULL);
			else
				continue;
			used_fields++;
		}

		// missing fields
		if (used_fields != number_fields)
			errx(7, "Impossible to read some of the data from file");
	}

	if (i != n_atoms)
		errx(8, "Loaded %ld/%ld atoms\n", i, n_atoms);

	fclose(handle);
}

void write_out(FILE *handle, simulation_space *space, int time, double kinetic, double potential)
{
	// using %g's for compatibility with lammps output
	fprintf(handle, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%ld\nITEM: BOX BOUNDS %s %s %s\n%g %g\n%g %g\n%g %g\nITEM: ATOMS id mass x y z vx vy vz fx fy fz\n",
	        time, space->number_atoms, space->period.x ? "pp" : "ff", space->period.y ? "pp" : "ff", space->period.z ? "pp" : "ff",
	        space->lo.x, space->hi.x, space->lo.y, space->hi.y, space->lo.z, space->hi.z);

	int i;
	for (i = 0; i < space->number_atoms; i++)
		fprintf(handle, "%ld %g %g %g %g %g %g %g %g %g %g\n", id[i], m[i], x[i].x, x[i].y, x[i].z, v[i].x, v[i].y, v[i].z, f[i].x, f[i].y, f[i].z);

	fprintf(handle, "ITEM: POTENTIAL ENERGY\n%.15e\nITEM: KINETIC ENERGY\n%.15e\n", potential, kinetic);
}

void write_file(char *filename, simulation_space *space, int time, double kinetic, double potential)
{
	FILE *handle = fopen(filename, "w");
	write_out(handle, space, time, kinetic, potential);
	fclose(handle);
}

