#!/bin/bash

# call from wrapper script
# last updated 08/22/2022 by M. Clark

echo "before calling source: $PATH"
source activate slim
echo "after calling source: $PATH"

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SLiM/3.7.1
module list

today_date=$(date +%m%d%y)
# define variables using array key
echo this is task ${SLURM_ARRAY_TASK_ID}

read -r AVG_AGE R REP N BURN <<< $(awk -v row=${SLURM_ARRAY_TASK_ID} -F',' 'NR==row {print $1, $2, $3, $4, $5}' "$ARRAY_KEY" | tr -d '\r')

#convert to numbers
AVG_AGE=$((AVG_AGE))
R=$((R))
REP=$((REP))
N=$((N))
BURN=$((BURN))


# print variables imported from array key
echo ----------------------------------------------------------------------------------------
echo variables imported from $ARRAY_KEY
echo AVG_AGE is $AVG_AGE
echo R is $R
echo REP is $REP
echo N is $N
echo BURN is $BURN 

# print vital variables imported from wrapper
echo ----------------------------------------------------------------------------------------
echo the date is $DATE but running on ${today_date}
echo runname is $RUN_NAME
echo jobname is $JOBNAME
echo slimscript is $SLIMSCRIPT
echo This file is $EXECUTABLE
echo Logs for this job are in $LOGFILESDIR
echo logs are called: ${JOBNAME}_${SIM_BLOCK}_%A_${SLURM_ARRAY_TASK_ID}.out and ${JOBNAME}_${SIM_BLOCK}_%A_${SLURM_ARRAY_TASK_ID}.err
echo
echo
echo datafiles created by this job will be moved to ${HOMEDIR}/${OUTDIR}
echo output will use the header ${JOBNAME}_${AVG_AGE}_${R}_${REP}
echo running with $CPUS
echo ----------------------------------------------------------------------------------------

# move everything to TMPDIR
echo navigating to $TMPDIR
cd $TMPDIR

if [ ! -d ./${OUTDIR} ]; then mkdir ./${OUTDIR}; fi


echo copying ${INDIR}/${SLIMSCRIPT} to $TMPDIR
cp ${INDIR}/${SLIMSCRIPT} $TMPDIR

echo I am running ${SLIMSCRIPT} and putting the output in ${OUTDIR} 

slim -d N=${N} -d AVG_AGE=${AVG_AGE} -d R=${R} -d Q=${BURN} -d "O='${OUTDIR}/tree_${JOBNAME}_${AVG_AGE}_${R}_${REP}.trees'" -d "MI='${OUTDIR}/metaInd_${JOBNAME}_${AVG_AGE}_${R}_${REP}.txt'" -d "MA='${OUTDIR}/metaAll_${JOBNAME}_${AVG_AGE}_${R}_${REP}.txt'" ${SLIMSCRIPT}

echo done running the simulation

cp -r $OUTDIR $HOMEDIR

echo done copying output back to $HOMEDIR

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
