#!/bin/bash

echo "before calling source: $PATH"
source activate slim
echo "after calling source: $PATH"

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SLiM/2021dev
module list

echo running with N == ${N}
echo running with prob of mortality == ${P}

slim -d P=${P} -d N=${N} -d "M='${OUTDIR}/${HEADER}_${JOBNAME}_${REP}.txt'" ${SLIMSCRIPT}

echo done running the simulation

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------


