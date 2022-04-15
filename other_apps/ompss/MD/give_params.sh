#!/bin/bash

if [[ $# -lt 2 || "${1%.lammps}" == "$1" || ("${2%.lammps}" == "$2" && "${2%.in}" == "$2") || ! -f $1 || ! -f $2 ]]; then
	echo "Gives paramters for a run."
	echo "Alternately, wraps the given executable with correct parameters, called with execve"
	echo "Usage : $(basename $0) params.lammps dumpfile.in [executable [out_file]]";
	echo "out_file defaults to ./results/dumpfile.out"
	exit 1
fi

lammps=$1
dump=$2
shift 2


cutoff=(`grep "lj/cut" $lammps`)
cutoff=${cutoff[2]}

iterations=(`grep "^run" $lammps`)
iterations=${iterations[1]}

timestep=(`grep "^timestep" $lammps`)
timestep=${timestep[1]}

write_freq=(`grep "^dump" $lammps`)
write_freq=${write_freq[4]}

pair_style=(`grep "^pair_coeff" $lammps | sed 's/\*/"*"/g'`)
epsilon="${pair_style[3]}"
sigma="${pair_style[4]}"

#n_atoms=$(n sed -n '/ITEM: NUMBER OF ATOMS/{N;s/.*\n//p;q}' $dump)

ARGS="-epsilon=$epsilon -sigma=$sigma -cutoff=$cutoff -it=$iterations -dt=$timestep -write_freq=$write_freq $dump"

if (( ! $# )); then
	echo $ARGS
	exit 0
else
	exe=$1
	shift
fi

if (( $# )); then
	ARGS+="-write_file=$1"
	shift
fi

exec $exe $ARGS
