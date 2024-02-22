#!/bin/bash
		
# Last updated 01/31/2024 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-run_processing [model] [avg_age]
# command line variables: 
model=$1 # nWF or pWF from command line
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

# tree processing variables: 
mu=1e-8

#if [[ $model == nWF ]]
#then
	#if [[ $avg_age == 2 ]]; then gen=3; fi
        #if [[ $avg_age == 5 ]]; then gen=6; fi 
        #if [[ $avg_age == 10 ]]; then gen=11; fi 
        #if [[ $avg_age == 20 ]]; then gen=21; fi 
#	if [[ $avg_age == 2 ]]; then gen=2.999165; fi
#	if [[ $avg_age == 5 ]]; then gen=6.006728; fi 
#	if [[ $avg_age == 10 ]]; then gen=11.01567; fi 
#	if [[ $avg_age == 20 ]]; then gen=20.97633; fi 
#else 
#	gen=1 # for pWF model
#fi

# define upper level variables:
jobname=run-trees #label for SLURM book-keeping 
run_name=DC_slim #label to use on output files

if [[ $model == nWF ]]
then
    if [[ $array_key == high ]]
    then
	date=02102024
    fi
    
    if [[ $array_key == low ]]
    then
	date=02112024
    fi
    treeprocess=tree_2_sum.py #processing python script 
fi

if [[ $model == pWF ]]
then
    date=04262023
    header=${model}
    treeprocess=pWF_tree_2_sum.py #processing python script 
fi

rundate=$(date +%m%d%Y)
# define dirs:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/py_logfiles_sumstats_${date} # name of directory to create and then write log files to
indir=$storagenode/$run_name/slim_output_${date} # where tree files live
pythondir=$storagenode/$run_name/scripts # where the python file lives
outdir=sum_stat_output_${rundate}
homedir=$storagenode/$run_name/

# define files
executable=$storagenode/$run_name/scripts/run_processing.sbatch #script to run 

# running variables
cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=8G #amount of RAM to request/use per CPU 
reps=100 # testing with 1 rep, should be 100 

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d $logfilesdir ]; then mkdir $logfilesdir; fi
if [ ! -d $outdir ]; then mkdir ./$outdir; fi

#submit job to cluster
				
sbatch --job-name=$jobname \
	--array=1-600 \
	--export=JOBNAME=$jobname,TREEPROCESS=$treeprocess,ARRAY_KEY=$array_key,MODEL=$model,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,HOMEDIR=$homedir,PYTHONDIR=$pythondir,MU=$mu,DATE=$date,EXECUTABLE=$executable,HEADER=$header,REPS=$reps,LOGFILESDIR=$logfilesdir \
	--cpus-per-task=$cpus \
	--mem-per-cpu=$ram_per_cpu \
	--output=$logfilesdir/${jobname}_${sim_block}_${date}_%A_%A.out \
	--error=$logfilesdir/${jobname}_${sim_block}_${date}_%A_%a.err \
	--time=24:00:00 \
	$executable

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	
