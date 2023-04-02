#!/bin/bash
		
# Last updated 11/14/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

# usage: ./scripts/wrapper-run_processing [model] [avg_age]
# command line variables: 
model=$1 # nWF or pWF from command line
avg_age=$2

# tree processing variables: 
mu=1e-8

if [[ $model == nWF ]]
then
	if [[ $avg_age == 2 ]]; then gen=2.999165; fi
	if [[ $avg_age == 5 ]]; then gen=6.006728; fi 
	if [[ $avg_age == 10 ]]; then gen=11.01567; fi 
	if [[ $avg_age == 20 ]]; then gen=20.97633; fi 
else 
	gen=1.00 # for pWF model
fi

# define upper level variables:
jobname=run-trees #label for SLURM book-keeping 
run_name=DC_slim #label to use on output files

if [[ $model == nWF ]]
then
    date=11082022 
    header=${model}_${avg_age} # from input when running wrapper-run_slim_all.sh
fi

if [[ $model == pWF ]]
then
    date=05312022
    header=${model}
fi
# define dirs:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
logfilesdir=$storagenode/$run_name/py_logfiles_sumstats_${date} # name of directory to create and then write log files to
indir=$storagenode/$run_name/slim_output_${date} # where tree files live
pythondir=$storagenode/$run_name/scripts # where the python file lives
outdir=sum_stat_output_${date}
homedir=$storagenode/$run_name/

# define files
executable=$storagenode/$run_name/scripts/run_processing.sbatch #script to run 
treeprocess=tree_2_sum.py #processing python script 

# running variables
cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=4G #amount of RAM to request/use per CPU 
reps=100 # testing with 1 rep, should be 100 

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d $logfilesdir ]; then mkdir $logfilesdir; fi
if [ ! -d $outdir ]; then mkdir ./$outdir; fi

    #submit job to cluster
for r in 2 10 100; do
	for rep in $(seq 1 $reps) ; do 
        if [[ $model == nWF ]]
        then
             filename=tree_${model}_${avg_age}_${r}_${rep}.trees
             metafile=metaInd_${model}_${avg_age}_${r}_${rep}.txt
             output_file=${homedir}${outdir}/${model}_${avg_age}_${r}_${rep}_${date}_age_bins.txt
        fi
        
        if [[ $model == pWF ]]
		then
             filename=tree_${model}_${r}_${rep}.trees
             metafile=metaInd_${model}_${r}_${rep}.txt
             output_file=${homedir}${outdir}/${model}_${r}_${rep}_${date}_age_bins.txt
        fi

		if [ ! -f "$output_file" ] # don't start job if output files already exists
		then 
				
			echo Starting job ${model}_${r}_${rep}_${date} with $filename and $metafile
				
			sbatch --job-name=$jobname \
			--export=JOBNAME=$jobname,TREEPROCESS=$treeprocess,MODEL=$model,FILENAME=$filename,METAFILE=$metafile,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,HOMEDIR=$homedir,PYTHONDIR=$pythondir,MU=$mu,R=$r,AVG_AGE=$avg_age,GEN=$gen,DATE=$date,EXECUTABLE=$executable,HEADER=$header,REPS=$reps,LOGFILESDIR=$logfilesdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=$logfilesdir/${header}_${r}_${rep}_${date}_%A.out \
			--error=$logfilesdir/${header}_${r}_${rep}_${date}_%A.err \
			--time=12:00:00 \
			$executable
		else
			echo output files for ${model}_${r}_${rep}_${date} already exist in ${outdir}
		fi		
	done
done	

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	
