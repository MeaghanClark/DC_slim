#!/bin/bash --login

########## SBATCH Lines for Resource Request ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --cpus-per-task=1         # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=12G            # memory required per allocated CPU (or core)
#SBATCH --job-name=slim_single      # you can give your job a name for easier identification (same as -J)
#SBATCH --output="/mnt/home/clarkm89/DC_slim/slim_log_single_07272023.out"
#SBATCH --error="/mnt/home/clarkm89/DC_slim/slim_log_single_07272023.err"

##########

# call from wrapper script
# last updated 07/27/2023 by M. Clark

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

# define variables
HOMEDIR=/mnt/home/clarkm89/DC_slim
OUTDIR=slim_output_troubleshooting_07272023
INDIR=$HOMEDIR/slim
SLIMSCRIPT="troubleshooting_nWF.slim"

# move everything to TMPDIR
echo navigating to $TMPDIR
cd $TMPDIR

if [ ! -d ./${OUTDIR} ]; then mkdir ./${OUTDIR}; fi

echo copying ${INDIR}/${SLIMSCRIPT} to $TMPDIR
cp ${INDIR}/${SLIMSCRIPT} $TMPDIR


echo I am running ${SLIMSCRIPT} and putting the output in ${OUTDIR} 

slim -d "O='${OUTDIR}/nWF_troubleshooting.trees'" ${SLIMSCRIPT}

echo done running the simulation

cp -r $OUTDIR $HOMEDIR

echo done copying output back to $HOMEDIR

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

echo ----------------------------------------------------------------------------------------
