#!/bin/bash
		
# Last updated 07/28/2021 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)
# usage: ./scripts/wrapper-run_slim_all.sh [pWF or nWF]

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

# ***JOB NAME IS IN COMMAND LINE ARGUMENT*** 
jobname=$1 #label for SLURM book-keeping, nWF or pWF 
run_name=DC_slim #label to use on output files
logfilesdir=logfiles #name of directory to create and then write log files to
executable=$storagenode/$run_name/scripts/run_slim_all.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=$storagenode/$run_name/slim/demo_change_${jobname}.slim #slimulation to run
treeprocess=$storagenode/$run_name/scripts/DC_slim_hpcc.py #processing python script

kvalue=7500
bvalue=100
outdir=$storagenode/$run_name/output

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=75G #amount of RAM to request/use per CPU 



#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
for rep in {1..10} ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,TREEPROCESS=$treeprocess,KVALUE=$kvalue,BVALUE=$bvalue,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
		--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
		--time=48:00:00 \
		$executable
done	

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

