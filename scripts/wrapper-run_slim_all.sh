#!/bin/bash
		
# Last updated 04/05/2023 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-run_slim_het_test.sh [pWF or nWF] [Nc] [avg_age] [R]

# command line variables: 
jobname=$1 #label for SLURM book-keeping, nWF or pWF 

# slim specific variables
n=$2 # census pop size
reps=1 # reps of slimulation to run 
avg_age=$3
r=$4 # factor to reduce pop size by

# select burn-in based on avg_age (Ne will change with age, so burn in must also change!) burn-in = avg_age * 10 * Ne
	# Ne estimated in estimate_Ne.R

# 07292023: trying avg_Age * 2 * Ne instead
if [[ $avg_age == 2 ]]
then
#	burn=126260
#	burn=126300
	burn=25300
elif [[ $avg_age == 5 ]]
then
#	burn=315890
#	burn=315900
	burn=63200
elif [[ $avg_age == 10 ]] 
then
#	burn=635576
#	burn=635600
	burn=127150
elif [[ $avg_age == 20 ]] 
then
#	burn=1264558
#	burn=1264600
	burn=253000

elif [[ $avg_age == 1 ]]
then
	burn=100000
else
	echo average age is invalid
fi


#define upper-level variables:
date=$(date +%m%d%Y)
#date=04062023
header=${jobname}_${avg_age}_${r} # header name, can change
run_name=DC_slim #label to use on output files

# define dirs
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/slim_log_${date} #name of directory to create and then write log files to
outdir=slim_output_${date}
indir=$storagenode/$run_name/slim
homedir=$storagenode/$run_name

executable=$storagenode/$run_name/scripts/run_slim_all.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=demo_change_${jobname}.slim #slimulation to run

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=12G #amount of RAM to request/use per CPU 

if [[ $avg_age == 2 ]]
then
	time=12:00:00
elif [[ $avg_age == 5 ]]
then
	time=24:00:00
elif [[ $avg_age == 10 ]] 
then
	time=48:00:00
elif [[ $avg_age == 20 ]] 
then
	time=168:00:00
elif [[ $avg_age == 1 ]] 
then
	time=12:00:00
else
	echo average age is invalid
fi

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d $logfilesdir ]; then mkdir $logfilesdir; fi

#submit job to cluster
	for rep in $(seq 1 $reps) ; do 
#	for rep in 61 ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,DATE=$date,SLIMSCRIPT=$slimscript,N=$n,AVG_AGE=$avg_age,BURN=$burn,P=$p,R=$r,HEADER=$header,REPS=${reps},REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,INDIR=$indir,HOMEDIR=$homedir,LOGFILESDIR=$logfilesdir,EXECUTABLE=$executable \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=$logfilesdir/${header}_${rep}_%A.out \
		--error=$logfilesdir/${header}_${rep}_%A.err \
		--time=$time \
		$executable
		
		echo submitting job with an average age of $avg_age and N of $n to run for $time!
	done	
#done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

