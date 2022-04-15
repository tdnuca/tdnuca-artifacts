# Molecular Dynamics @ OmpSs

A simple Molecular Dynamics (MD) simulation using the Lennard-Jones potential.

Author: Luc Jaulmes <luc.jaulmes@bsc.es><br />
Application written during summer school at LLNL, released under open source (details of license TBD).

---

This code implements a simple MD benchmark, simulating atoms interacting based on a Lennard-Jones potential.

The implementation is Arrays-of-Structs style. The code subdivides the simulation space into boxes based on the cutoff
distance of the force and assigns atoms to boxes by sorting the global arrays of atoms and assigning slices of this
array to each box.

The program takes as input a dump file from LAMMPS listing all the atoms and some supplementary parameters on the
command line (which can be extracted from a LAMMPS input file by using the `give_params.sh` script). Its output is a
file that outputs the atoms data in the same format, with the addition of total kinetic and potential energies in the
system.

Features:
- Data-flow dependencies
- `reduction:+` for force computations, that can be changed to `commutative`, or to `concurrent` provided the udpates to
  `f` and `pe` in `single_force()` are made atomic

What *really needs* to be improved:
- Allow tasks to be at a coarser granularity than single boxes.

What *could* be improved:
- Skip energy calculations on non-output iterations, in particular implement leapfrog integration for velocities
- Parallelize sorting of atoms
- Increase box size to increase the number of iterations between atom sorts
  (This implies considering more interactions that are beyond the cutoff distance than with smalelst box size!)
- Output energies withotu dumping all the atoms

What can be changed:
- Each box could be a separate array instead of a slice of the global array
	* (+) constant dependency addresses from one iteration to the next
	* (-) keep track of how much space is allocated per box, make sure it is always enough
	* (-) complicates sorting


## Get started
Be sure to have OmpSs set up in your path and library path.

To compile, simply cd to directory and type make. This will produce 4 executables
- `md_perf` : your parallel MD simulation
- `md_inst` : same program linked with the instrumentation version of the OmpSs runtime, allows tracing etc.
- `md_debug`: same program linked with the debug version of the OmpSs runtime, and compiled with -O0
- `md_seq`  : sequential version of the MD simulation (compiled by ignoring pragmas)

Then run with `./md_perf [options] lammps.dump`, or easier:

    ./md_perf $(./give_params.sh lammps.in lammps.dump)

Parameters are read from the command line, and are:

| Parameter       | in LAMMPS configuration            |   |
| --------------- | ---------------------------------- | - |
|`-dt=0.005`      | `timestep <value>`                 | The timestep, i.e. the duration of an iteration <br />(thus frequency at which we compute forces, positions, velocities) |
|`-epsilon=1.0`   | `pair_coeff * * <epsilon> <sigma>` | Parameter of the LJ potential |
|`-sigma=1.0`     | `pair_coeff * * <epsilon> <sigma>` | Parameter of the LJ potential |
|`-cutoff=2.5`    | `pair_style lj/cut <value>`        | Distance between particles beyond which interactions are not computed |
|`-it=100`        | `run <value>`                      | The number of iterations to run |
|`-write_freq=100`| `dump <options>`                   | Frequency (in iterations) to dump outputs and energies |
|`-write-file=`   | `dump <options>`                   | Output file for atoms dump, defaults to stdout |

Parameters can be in any unit, as long as the system is consistent, since the program doesn't make any conversions, for
example all LJ units. If units used are Length in `Angstrom`, Time in `ps`, Energy in `eV`, then masses need to be
`eV*ps^2/Ang^2`.

## Inputs

Inputs are available in `inputs.tgz` and explained in the contained `README_inputs.tgz`

## Additional information
You can use a visualization tool such as [Ovito](https://ovito.org/) to see dump files and animate what happens over time.

The `give_params.sh` can also wrap the executable which is useful for e.g. gdb:

    $ gdb ./md_debug -ex "set exec-wrapper ./give_params.sh lammps.in lammps.dump"
	# to run printing atoms to stdout
	(gdb) run
	# to runs printing atoms to /dev/null
	(gdb) run /dev/null
