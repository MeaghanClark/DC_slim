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
elif [[ $sim_block == "WF" ]]
then
	array_key=/mnt/home/clarkm89/DC_slim/scripts/array_index_key_1.txt
else
	echo simulation type invalid
fi
# slim specific variables
#n=61525 # census op size
#res=100 # reps of slimulation to run 
#avg_age=20
#r=100 # factor to reduce op size by

# select burn-in based on avg_age (Ne will change with age, so burn in must also change!) burn-in = avg_age * 10 * Ne
	# Ne estimated in estimate_Ne.R

# 07292023: trying avg_Age * 2 * Ne instead
# 01282024: will secify in array key
#if [[ $avg_age == 2 ]]
#then
#	burn=126260
#	burn=126300
#	burn=25300
#	burn=120000 # 01282024, for Ne of around 30000
#elif [[ $avg_age == 5 ]]
#then
#	burn=315890
#	burn=315900
#	burn=63200
#	burn=300000
#elif [[ $avg_age == 10 ]] 
#then
#	burn=635576
#	burn=635600
#	burn=127150
#	burn=600000
#elif [[ $avg_age == 20 ]] 
#then
#	burn=1264558
#	burn=1264600
#	burn=253000
#	burn=1200000

#elif [[ $avg_age == 1 ]]
#then
#	burn=100000
#else
#	echo average age is invalid
#fi


#define uper-level variables:
date=$(date +%m%d%Y)
# header=${jobname}_${avg_age}_${r} # header name, can change
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

#submit job to cluster
	# for re in $(seq 1 $reps) ; do 
		sbatch --job-name=$jobname \
		--array=1-600 \
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
	

