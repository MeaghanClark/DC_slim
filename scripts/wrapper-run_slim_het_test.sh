#!/bin/bash
		
# Last updated 04/21/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)
# usage: ./scripts/wrapper-run_slim_het_test.sh [pWF or nWF]

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

# ***JOB NAME IS IN COMMAND LINE ARGUMENT*** 
jobname=$1 #label for SLURM book-keeping, nWF or pWF 
run_name=DC_slim #label to use on output files
logfilesdir=logfiles_het_test #name of directory to create and then write log files to
executable=$storagenode/$run_name/scripts/run_slim_het_test.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=$storagenode/$run_name/slim/het_test_demo_change_${jobname}.slim #slimulation to run
outdir=$storagenode/$run_name/het_test_output

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=12G #amount of RAM to request/use per CPU 

n=7500
p=0.2

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
	for rep in {1..100} ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,N=$n,P=$p,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${jobname}_${rep}_%A.out \
		--error=./$logfilesdir/${jobname}_${rep}_%A.err \
		--time=48:00:00 \
		$executable
		
		echo submitting job with prob of mortality of $p and N of $n!
	done	
#done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

