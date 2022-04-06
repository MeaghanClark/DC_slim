#!/bin/bash

#--------------- EXECUTABLE ---------------

# 		--export=JOBNAME=$jobname,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \

source activate slim


#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module list 


# sys.argv = ['tree_processing.py', './path/to/tree/file', './outdir', 'big_test']

python -u ./scripts/Get_ped_ids.py ./output/nWF_${REP}/tree_nWF_${REP}.trees ./${OUTDIR}/nWF_${REP}/labels/ nWF_${REP}

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------

