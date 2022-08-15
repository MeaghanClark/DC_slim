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

today_date=$(date +%m%d%y)

# print vital variables imported from wrapper
echo ----------------------------------------------------------------------------------------
echo the date is $DATE but running on ${today_date}
echo runname is $RUN_NAME
echo jobname is $JOBNAME
echo slimscript is $SLIMSCRIPT
echo This file is $EXECUTABLE
echo Logs for this job are in $LOGFILESDIR
echo logs are called: ${HEADER}_${REP}_%A.out and ${HEADER}_${REP}_%A.err
echo
echo
echo This is rep ${REP} / ${REPS} of ${HEADER}
echo running with N == ${N}
echo running with AVG_AGE == ${AVG_AGE}
echo running with R == ${R}
echo
echo
echo datafiles created by this job will be moved to ${HOMEDIR}/${OUTDIR}
echo output will use the header $HEADER
echo running with $CPUS
echo ----------------------------------------------------------------------------------------

# move everything to TMPDIR
echo navigating to $TMPDIR
cd $TMPDIR

if [ ! -d ./${OUTDIR} ]; then mkdir ./${OUTDIR}; fi


echo copying ${INDIR}/${SLIMSCRIPT} to $TMPDIR
cp ${INDIR}/${SLIMSCRIPT} $TMPDIR

echo I am running ${SLIMSCRIPT} and putting the output in ${OUTDIR} 

slim -d AVG_AGE=${AVG_AGE} -d N=${N} -d "M='${OUTDIR}/${HEADER}_${REP}.txt'" ${SLIMSCRIPT}

cp -r $OUTDIR $HOMEDIR

echo done copying output back to $HOMEDIR

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------


