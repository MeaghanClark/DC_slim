#!/bin/bash
		
# Last updated 04/28/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

jobname=run-trees_inc_mu_test #label for SLURM book-keeping 
run_name=DC_slim #label to use on output files
logfilesdir=logfiles_trees_to_pi #name of directory to create and then write log files to
executable=$storagenode/$run_name/scripts/run_processing.sbatch #script to run 
treeprocess=$storagenode/$run_name/scripts/tree_2_pi_inc_mu.py #processing python script

outdir=$storagenode/$run_name/het_inc_mu
indir=$storagenode/$run_name/full_output

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=500G #amount of RAM to request/use per CPU 

model=nWF # nWF or pWF from command line
header=inc_mu_52322 # from input when running wrapper-run_slim_all.sh

#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

    #submit job to cluster
#    for rep in {1..10} ; do 
        filename=${indir}/tree_${model}_full_run_42922_74.trees
	rep=74
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,TREEPROCESS=$treeprocess,MODEL=$model,FILENAME=$filename,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${jobname}_${model}_${rep}_%A.out \
		--error=./$logfilesdir/${jobname}_${model}_${rep}_%A.err \
		--time=168:00:00 \
		$executable
#    done	

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

