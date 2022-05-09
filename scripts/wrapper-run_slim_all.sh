#!/bin/bash
		
# Last updated 04/27/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)
# usage: ./scripts/wrapper-run_slim_het_test.sh [pWF or nWF] [Nc] [header]

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

# ***JOB NAME IS IN COMMAND LINE ARGUMENT*** 
jobname=$1 #label for SLURM book-keeping, nWF or pWF 
run_name=DC_slim #label to use on output files
logfilesdir=logfiles #name of directory to create and then write log files to
executable=$storagenode/$run_name/scripts/run_slim_all.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=$storagenode/$run_name/slim/demo_change_${jobname}.slim #slimulation to run
outdir=$storagenode/$run_name/full_output_43022

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=75G #amount of RAM to request/use per CPU 

n=$2
p=0.2
header=$3

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
	for rep in {1..10} ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,N=$n,P=$p,HEADER=$header,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${header}_${rep}_%A.out \
		--error=./$logfilesdir/${header}_${rep}_%A.err \
		--time=24:00:00 \
		$executable
		
		echo submitting job with prob of mortality of $p and N of $n!
	done	
#done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

