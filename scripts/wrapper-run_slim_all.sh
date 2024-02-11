#!/bin/bash
		
# Last udated 01/29/2024 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scrits/wrapper-run_slim_het_test.sh [pWF or nWF] [Nc] [avg_age] [R]

# command line variables: 
jobname=$1 #label for SLURM book-keeing, nWF or pWF 
sim_block=$2
if [[ $sim_block == "low" ]]
then
	array_key=/mnt/home/clarkm89/DC_slim/scripts/array_index_key_2_5.txt
elif [[ $sim_block == "high" ]]
then
	array_key=/mnt/home/clarkm89/DC_slim/scripts/array_index_key_10_20.txt
elif [[ $sim_block == "pWF" ]]
then
	array_key=/mnt/home/clarkm89/DC_slim/scripts/array_index_key_1.txt
else
	echo simulation type invalid
fi

#define uper-level variables:
date=$(date +%m%d%Y)
run_name=DC_slim #label to use on outut files

# define dirs
storagenode=/mnt/home/clarkm89 #ath to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/slim_log_${date} #name of directory to create and then write log files to
outdir=slim_output_${date}
indir=$storagenode/$run_name/slim
homedir=$storagenode/$run_name

executable=$storagenode/$run_name/scripts/run_slim_sum.sbatch #script to run 


# *** CHANGE SLIM SCRIPT BASED ON COMMAND LINE***
slimscript=demo_change_${jobname}.slim #slimulation to run

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=8G #amount of RAM to request/use per CPU 

if [[ $sim_block  == "low" ]]
then
	time=72:00:00
elif [[ $sim_block == "high" ]]
then
	time=72:00:00
elif [[ $sim_block == "pWF" ]] 
then
	time=24:00:00
else
	echo sim_block is invalid
fi

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d $logfilesdir ]; then mkdir $logfilesdir; fi


job_no=$(wc -l < "$array_key")
echo "The number of jobs we are requestion is $job_no"

#submit job to cluster
	# for re in $(seq 1 $reps) ; do 
		sbatch --job-name=$jobname \
		--array=1-$job_no \
		--export=JOBNAME=$jobname,DATE=$date,SLIMSCRIPT=$slimscript,ARRAY_KEY=$array_key,SIM_BLOCK=$sim_block,HEADER=$header,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,INDIR=$indir,HOMEDIR=$homedir,LOGFILESDIR=$logfilesdir,EXECUTABLE=$executable \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=$logfilesdir/${jobname}_${sim_block}_%A_%a.out \
		--error=$logfilesdir/${jobname}_${sim_block}_%A_%a.err \
		--time=$time \
		$executable
		
		echo submitting job using $array_key! 
	#done	
#done

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

