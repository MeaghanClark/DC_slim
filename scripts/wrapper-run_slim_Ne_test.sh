#!/bin/bash
		
# Last updated 09/26/2022 by MI Clark, script format by R Toczydlowski 

# run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-run_slim_Ne_test.sh [N] [avg_age]

# command line variables: 

# slim specific variables 
jobname=Ne_test
n=$1 # census pop size
reps=100 # reps of slimulation to run 
avg_age=$2

# define upper-level variables:
date=$(date +%m%d%Y)
header=${jobname}_${avg_age}_${n} # header name, can change
run_name=DC_slim #label to use on output files

# define dirs
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/slim_logfiles_Ne_test_${date} #name of directory to create and then write log files to
outdir=Ne_test_output_${date}
indir=$storagenode/$run_name/slim
homedir=$storagenode/$run_name

executable=$storagenode/$run_name/scripts/run_slim_Ne_test.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=$storagenode/$run_name/slim/Ne_testing_nWF.slim #slimulation to run

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=2G #amount of RAM to request/use per CPU 

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d $logfilesdir ]; then mkdir $logfilesdir; fi

#submit job to cluster
	for rep in $(seq 1 $reps) ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,HEADER=$header,N=$n,AVG_AGE=$avg_age,REP=$rep,REPS=$reps,CPUS=$cpus,RUN_NAME=$run_name,DATE=$date,STORAGENODE=$storagenode,OUTDIR=$outdir,HOMEDIR=$homedir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=$logfilesdir/${jobname}_${avg_age}_${n}_${rep}_%A.out \
		--error=$logfilesdir/${jobname}_${avg_age}_${n}_${rep}_%A.err \
		--time=4:00:00 \
		$executable
		
		echo submitting job with an average age of $avg_age and N of $n!
	done	
#done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

