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

today_date=$(date +%m%d%y)

# print vital variables imported from wrapper
echo ----------------------------------------------------------------------------------------
echo the date is $DATE but running on ${today_date}
echo runname is $RUN_NAME
echo jobname is $JOBNAME
echo This file is $EXECUTABLE
echo Logs for this job are in $LOGFILESDIR
echo logs are called: ${HEADER}_${R}_${REP}_%A.out and ${HEADER}_${R}_${REP}_%A.err
echo
echo
echo This is rep ${REP} / ${REPS} of ${HEADER}
echo running with mu == ${MU}
echo running with gen == ${GEN}
echo
echo
echo datafiles created by this job will be moved to ${HOMEDIR}/${OUTDIR}
echo output will use the prefix ${MODEL}_${REP}
echo running with $CPUS
echo ----------------------------------------------------------------------------------------


# move everything to tmpdir
cd $TMPDIR
if [ ! -d ./$OUTDIR ]; then mkdir ./$OUTDIR; fi

# move python script
cp ${PYTHONDIR}/${TREEPROCESS} $TMPDIR

# move tree file
cp ${INDIR}/${FILENAME} $TMPDIR

# arguments: 
# [0] -- python script name
# [1] -- tree file
# [2] -- outdir
# [3] -- prefix
# [4] -- mu
# [5] -- gen time

# python -u ./pythonScript treefile outdir prefix mu gen_time
python -u ${TREEPROCESS} ${FILENAME} ${OUTDIR} ${MODEL}_${R}_${REP}_${DATE} ${MU} ${GEN} > ${LOGFILESDIR}/log_${JOBNAME}_${MODEL}_${R}_${REP}.out

echo done processing trees 

cp -r $OUTDIR $HOMEDIR

echo done copying files over to homedir

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------


