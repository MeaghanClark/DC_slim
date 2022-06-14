#!/bin/bash
		
# Last updated 04/28/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# command line variables: 
model=$1 # nWF or pWF from command line

# tree processing variables: 
mu=1e-8

if [[ $model == nWF ]]
then
	gen=5.00363475332601
else
	gen=1.00
fi

# define upper level variables:
jobname=run-trees #label for SLURM book-keeping 
run_name=DC_slim #label to use on output files
date=05312022 #$(date +%m%d%Y)
header=${model}_${date} # from input when running wrapper-run_slim_all.sh

# define dirs:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/py_logfiles_${date} #name of directory to create and then write log files to
indir=$storagenode/$run_name/slim_output_05312022 # where tree files live
pythondir=$storagenode/$run_name/scripts # where the python file lives
outdir=het_output_${date}
homedir=$storagenode/$run_name/

# define files
executable=$storagenode/$run_name/scripts/run_processing.sbatch #script to run 
treeprocess=tree_2_het.py #processing python script
# dataprefix=test

# running variables
cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=24G #amount of RAM to request/use per CPU 
reps=100

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d $logfilesdir ]; then mkdir $logfilesdir; fi

    #submit job to cluster
for r in 2 10 100; do  

	for rep in $(seq 1 $reps) ; do 
        	filename=tree_${model}_${r}_${rep}.trees
			sbatch --job-name=$jobname \
			--export=JOBNAME=$jobname,TREEPROCESS=$treeprocess,MODEL=$model,FILENAME=$filename,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,HOMEDIR=$homedir,PYTHONDIR=$pythondir,MU=$mu,R=$r,GEN=$gen,DATE=$date,EXECUTABLE=$executable,HEADER=$header,REPS=$reps,LOGFILESDIR=$logfilesdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=$logfilesdir/${header}_${r}_${rep}_%A.out \
			--error=$logfilesdir/${header}_${r}_${rep}_%A.err \
			--time=72:00:00 \
			$executable
	done
done	

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

