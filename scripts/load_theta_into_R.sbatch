#!/bin/bash


#load programs we want to use
# CHANGE
module purge
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load  R/4.0.2

module list

export R_LIBS_USER=~/R_Lib/4.0.2

echo I am importing: 

# run R script

# define args
# args[1] is indir, "../theta_output_dist_11082022"
# args[2] is outdir for Robj

Rscript ${RSCRIPT} --args ${INDIR} ${OUTDIR}

echo done running 

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------