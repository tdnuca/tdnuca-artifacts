## MD input files
> **These files and explanations are for running MD, courtesy of Sam Reeve <sreeve@purdue.edu>**

Most files are currently in LJ units(everything scaled to epsilon and sigma, easiest when epsilon = sigma = 1).
This includes everything from simulation box sizes, cutoff distance, timestep, positions, velocities, etc.
See http://lammps.sandia.gov/doc/units.html starting at "For style lj..."

Units whose filename contains "realunits" are consistent within LAMMPS except for mass(units style metal) as follows:
Length - Angstrom; Time - ps; Energy - eV; Mass - amu
Masses have been converted to eV\*ps^2/Ang^2 for a fully consistent set of units in all input files in this directory
(using `convert_realunits_mass.py`).

Files ending in ".lammps" are the input files for LAMMPS
Files ending in ".in" are what you should use to initialize your code.

Make sure the following parameters for LAMMPS values match with those you give your program:
- Interatomic potential is given by `pairstyle ****` (currently all lj/cut)
- Cutoff distance is given by `pair_style lj/cut ####` (default 2.5)
- Epsilon and sigma (in that order) are given by `pair_coeff * * #### ####` (default 1 and 1)
- The halo distance is given by `skin ####` (default 0.3)
- The duration of a timestep is given by `timestep ####` (default 0.001)
- The total number of iterations is given by `run ####` (default 1000)


Don't bother trying to match the following parameters, they're only necessary to produce the lists of atoms in LAMMPS:
  units, atom\_style, dimention, lattice, region, create\_box, mass, create\_atoms, velocity, dump, fix, unfix


### List of available problems
This file explains test simulations made for MD codes in various programming models. These examples are broken into
general categories with brief description of intention. All files are Lennard-Jones materials(epsilon and sigma = 1)
with Lennard-Jones reduced units unless otherwise specified as copper with real units(Angstroms, ps, eV, etc.)



| Toy problems               | |
| -------------------------- |-|
| 2atoms\_fff\_dump1         | atoms oscillate within energy basin |
| 2atoms\_ppp\_v5\_dump1     | atoms repel and interact around periodic boundary |
| 7atoms\_fff\_v3\_dump1     | forces perfectly cancel in each coordinate direction; low velocity so atoms stay inside simulation box |
| 2D\_5atoms\_fff\_v3\_dump1 | same as 7atoms\_fff\_v3\_dump1 in 2D |
| 7atoms\_ppp\_v5\_dump1     | forces perfectly cancel in each coordinate direction, both directly and around periodic boundaries |
| 2D\_5atoms\_ppp\_v5\_dump1 | same as 7atoms\_ppp\_v5\_dump1 in 2D |


| Frozen (0K)                      | |
| -------------------------------- |-|
| 1000atoms\_ppp\_dump10           | initial perfect simple cubic: unstable structure at zero Kelvin; round-off errors quickly break symmetry and atoms become disordered |
| 864atoms\_ppp\_dump10\_realunits | initial perfect FCC: stable structure at zero Kelvin; stays frozen for all simulations currently run |
| 2D\_800atoms\_ppf\_dump10        | same as 1000atoms\_ppp\_dump10 in 2D |


| LJ fluid                      | |
| ----------------------------- |-|
| 64atoms\_ppp\_v1\_dump1       | Small example of liquid with initial reduced temperature of 1 |
| 1000atoms\_ppp\_v1\_dump10    | Liquid with initial reduced temperature of 1 |
| 1000atoms\_ppp\_v3\_dump100   | Same as 1000atoms\_ppp\_v1\_dump10 with higher initial temperature |
| 2D\_800atoms\_ppf\_v1\_dump10 | Same as 1000atoms\_ppp\_v1\_dump10 in 2D |


| Copper                                   | |
| ---------------------------------------- |-|
| 864atoms\_ppp\_v300\_dump100\_realunits  | Copper FCC crystal with initial temperature of 300K |
| 864atoms\_ppp\_v4000\_dump100\_realunits | Same as 864atoms\_ppp\_v300\_dump100\_realunits, higher initial temp |
| 2048atoms\_ppp\_v300\_dump100\_realunits | Same as 864atoms\_ppp\_v300\_dump100\_realunits, larger |


| Nanoclusters               | |
| -------------------------- |-|
| 13atoms\_fff\_v1\_dump10   | Single cluster of 13 atoms, initially FCC |
| 55atoms\_fff\_v1\_dump10   | Single cluster of 55 atoms, initially FCC |
| 64atoms\_fff\_v1\_dump1    | Single cluster of 64 atoms, initially simple cubic (unstable) |
| 864atoms\_fff\_v1\_dump10  | Single cluster of 864 atoms, initially FCC |
| 351atoms\_fff\_v1\_dump100 | Multiple 13 atom clusters interacting |
| 440atoms\_fff\_v1\_dump100 | Multiple 55 atom clusters interacting |


| Phase change                         | |
| ------------------------------------ |-|
| 2000atoms\_ppp\_rho1.0\_T3\_dump100  | Reduced temperature of 3 and reduced density of 1; initially BCC structure changes phase to liquid (see LJ phase diagram) |
| 2000atoms\_ppp\_rho1.2\_T3\_dump100  | Reduced temperature of 3 and reduced density of 1.2; primarily keeps initial BCC structure |
| 2000atoms\_ppp\_rho1.4\_T3\_dump100  | Reduced temperature of 3 and reduced density of 1.4; changes phase from BCC to FCC |
| 2000atoms\_ppp\_rho1.6\_T3\_dump100  | Reduced temperature of 3 and reduced density of 1.6; changes phase from BCC to FCC; phase change faster than density 1.4 |
| 2000atoms\_ppp\_stretch\_T3\_dump100 | Reduced temperature of 3 and reduced density of 1 strained by 1% in one direction; similar phase change from BCC to FCC |
