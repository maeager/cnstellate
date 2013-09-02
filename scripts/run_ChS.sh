#!/bin/bash
#PBS -N run_ChS
#PBS -l walltime=24:0:0
#PBS -l pmem=8000MB
#PBS -A Monash038
# PBS -m ae


cd $PBS_O_WORKDIR
# cd `dirname $0`/..
module load octave/3.6.2

./mpi/special TStellate2.hoc -c "optimise_CS()"
# ./mpi/special AMResponses.hoc -c "SimpleResponses_CS()" -c "AMresponse_CS()"
# ./mpi/special AMResponses.hoc -c "AMresponse_CS()"

cd TStellate2_CS
../scripts/fullbattery.sh
