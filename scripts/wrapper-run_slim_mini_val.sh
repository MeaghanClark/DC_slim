#!/bin/bash
		
# Last updated 07/28/2021 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)
# usage: ./scripts/wrapper-run_slim_all.sh [pWF or nWF]

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

# ***JOB NAME IS IN COMMAND LINE ARGUMENT*** 
jobname=mini_val #label for SLURM book-keeping, nWF or pWF 
run_name=DC_slim #label to use on output files
logfilesdir=logfiles_minival #name of directory to create and then write log files to
executable=$storagenode/$run_name/scripts/run_slim_mini_val.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=$storagenode/$run_name/slim/mini_demo_change_nWF.slim #slimulation to run

# kvalues={2500 3500 4500 5500 6500 7500 8500 9500}
# bvalues are K/75
# fitness={0.5 0.6 0.7 0.8 0.9}

outdir=$storagenode/$run_name/mini_val_output

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=25G #amount of RAM to request/use per CPU 



#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
for kval in 2500 3500 4500 5500 6500 7500 8500 9500 ; do 
	for fitness_val in 0.5 0.6 0.7 0.8 0.9 ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,KVALUE=$kval,FITNESS=$fitness_val,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${jobname}_${kval}_${fitness_val}_%A.out \
		--error=./$logfilesdir/${jobname}_${kval}_${fitness_val}_%A.err \
		--time=12:00:00 \
		$executable
		
		echo submitting job with fitness scaling of $fitness_val and K-value of ${kval}!
	done	
done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

