#!/bin/bash

# usage: bash run_slim.sh file.slim K-value output

# importing: JOBNAME,SLIMSCRIPT,TREEPROCESS,KVALUE,BVALUE,REP,CPUS,RUN_NAME,STORAGENODE,OUTDIR,LOGFILESDIR \
echo "before calling source: $PATH"
source activate slim
echo "after calling source: $PATH"

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SLiM/2021dev
module list

echo This is rep ${REP} of ${HEADER}
echo running with N == ${N}
echo running with P == ${P}

if [ ! -d ${OUTDIR} ]; then mkdir ${OUTDIR}; fi

slim -d N=${N} -d P=${P} -d "O='${OUTDIR}/tree_${HEADER}_${REP}.trees'" -d "MI='${OUTDIR}/metaInd_${JOBNAME}_${REP}.txt'" -d "MA='${OUTDIR}/metaAll_${JOBNAME}_${REP}.txt'" ${SLIMSCRIPT}

echo done running the simulation

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
