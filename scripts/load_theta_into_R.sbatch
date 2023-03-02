#!/bin/bash


#load programs we want to use
# CHANGE
module purge
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load  R/4.0.2

module list

export R_LIBS_USER=~/R_Lib/4.0.2


# move everything to tmpdir
cd $TMPDIR
if [ ! -d ./$OUTDIR ]; then mkdir ./$OUTDIR; fi # make outdir in temp directory 
if [ ! -d ./$INDIR ]; then mkdir ./$INDIR; fi # make indir in temp directory 


# move R script
cp ${HOMEDIR}/scripts/${RSCRIPT} $TMPDIR
echo moved ${HOMEDIR}/scripts/${RSCRIPT} to $TMPDIR

# move tree file
cp ${HOMEDIR}/${INDIR}/* $TMPDIR/${INDIR}
echo moved files in ${HOMEDIR}/${INDIR}/ to $TMPDIR/${INDIR}, verified below

ls $TMPDIR/${INDIR} | wc -l 

# run R script

Rscript ${RSCRIPT} --args ${INDIR} ${OUTDIR}

echo done running 

cp -r $OUTDIR $HOMEDIR

echo done copying files over to homedir

ls ${HOMEDIR}/${OUTDIR} | wc -l 



#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
