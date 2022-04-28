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
module list

echo This is rep ${REP}

#export=JOBNAME=$jobname,TREEPROCESS=$treeprocess,MODEL=$model,FILENAME=$filename,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \

# python -u ./pythonScript treefile outdir prefix
python -u ${TREEPROCESS} ${FILENAME} ${OUTDIR} ${MODEL}_${REP} > ${LOGFILESDIR}/log_${JOBNAME}_${MODEL}_${REP}.out

echo done processing trees 

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------


