#!/bin/bash
		
# Last updated 08/15/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-run_slim_no_bottleneck.sh [pWF or nWF] [N] [avg_age]

# command line variables: 
jobname=$1 #label for SLURM book-keeping, nWF or pWF 

# slim specific variables 
n=$2 # census pop size
reps=10 # reps of slimulation to run 
avg_age=$3
p=1/(avg_age+1) # probability of mortality

# define upper-level variables:
date=$(date +%m%d%Y)
header=${jobname}_no_bln_${avg_age} # header name, can change
run_name=DC_slim #label to use on output files

# define dirs
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/slim_logfiles_no_bottleneck_${date} #name of directory to create and then write log files to
outdir=no_bottleneck_output_${date}
indir=$storagenode/$run_name/slim
homedir=$storagenode/$run_name

executable=$storagenode/$run_name/scripts/run_slim_no_bottleneck.sbatch #script to run 

# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=$storagenode/$run_name/slim/no_bottleneck_${jobname}.slim #slimulation to run

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=12G #amount of RAM to request/use per CPU 

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
	for rep in $(seq 1 $reps) ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,HEADER=$header,N=$n,P=$p,AVG_AGE=$avg_age,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,DATE=$date,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${jobname}_${rep}_%A.out \
		--error=./$logfilesdir/${jobname}_${rep}_%A.err \
		--time=12:00:00 \
		$executable
		
		echo submitting job with prob of mortality of $p and N of $n!
	done	
#done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

