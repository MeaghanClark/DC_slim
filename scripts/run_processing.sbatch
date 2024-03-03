#!/bin/bash

# usage: bash run_slim.sh file.slim K-value output
# importing: JOBNAME,SLIMSCRIPT,TREEPROCESS,KVALUE,BVALUE,REP,CPUS,RUN_NAME,STORAGENODE,OUTDIR,LOGFILESDIR \
echo "before calling source: $PATH"
source activate slim3
echo "after calling source: $PATH"

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
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
echo This file is $EXECUTABLE
echo Logs for this job are in $LOGFILESDIR
echo logs are called: ${HEADER}_${R}_${REP}_%A.out and ${HEADER}_${R}_${REP}_%A.err
echo
echo
echo This is rep ${REP} / ${REPS} of ${HEADER}
echo running with mu == ${MU}
echo running with gen == ${GEN}
echo Average age is ${AVG_AGE}
echo
echo datafiles created by this job will be moved to ${HOMEDIR}/${OUTDIR}
echo output will use the prefix ${MODEL}_${REP}
echo running with $CPUS
echo ----------------------------------------------------------------------------------------

# define input files 
FILENAME=tree_${MODEL}_${AVG_AGE}_${R}_${REP}.trees
METAFILE=metaInd_${MODEL}_${AVG_AGE}_${R}_${REP}.txt
OUTPUT_FILE=${HOMEDIR}${OUTDIR}/${MODEL}_${AVG_AGE}_${R}_${REP}_${DATE}_summary.txt # use if re-doing for subset of jobs.. not written yet! 

if [ -e "$OUTPUT_FILE" ]; then
       echo "Error: output file '$OUTPUT_FILE' already exists. Script terminated." 
       exit 1
else
       echo "Output file does not exist. Continuing script... " 
fi

if [[ $MODEL == nWF ]]
then
       if [[ $AVG_AGE == 2 ]]; then GEN=2.999165; fi
       if [[ $AVG_AGE == 5 ]]; then GEN=6.006728; fi 
       if [[ $AVG_AGE == 10 ]]; then GEN=11.01567; fi 
       if [[ $AVG_AGE == 20 ]]; then GEN=20.97633; fi 
else 
       GEN=1 # for pWF model
fi

# move everything to tmpdir
cd $TMPDIR
if [ ! -d ./$OUTDIR ]; then mkdir ./$OUTDIR; fi

# move python script
cp ${PYTHONDIR}/${TREEPROCESS} $TMPDIR

# move tree file
cp ${INDIR}/${FILENAME} $TMPDIR

# move metafile
cp ${INDIR}/${METAFILE} $TMPDIR


# arguments: 
# [0] -- python script name
# [1] -- tree file
# [2] -- meta file
# [3] -- outdir
# [4] -- prefix
# [5] -- mu
# [6] -- gen time
# [7] -- avg age
# [8] -- burn in 

python -u ${TREEPROCESS} ${FILENAME} ${METAFILE} ${OUTDIR} ${MODEL}_${AVG_AGE}_${R}_${REP}_${DATE} ${MU} ${GEN} ${AVG_AGE} ${BURN} ${R} > ${LOGFILESDIR}/log_${JOBNAME}_${MODEL}_${AVG_AGE}_${R}_${REP}.out

echo done processing trees 

cp -r $OUTDIR $HOMEDIR

echo done copying files over to homedir

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------


