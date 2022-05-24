#!/bin/bash
		
# Last updated 08/18/2020 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

#define variables:
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live

jobname=slim_Nc_test_WF #label for SLURM book-keeping 
run_name=DC_slim #label to use on output files
logfilesdir=logfiles_test #name of directory to create and then write log files to
executable=$storagenode/$run_name/scripts/run_slim_hpcc.sbatch #script to run 
slimscript=$storagenode/$run_name/slim/demo_change_WF.slim #slimulation to run
treeprocess=$storagenode/$run_name/scripts/DC_slim_hpcc.py #processing python script
kvalue=5000
bvalue=100
outdir=$storagenode/$run_name/output

cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=36G #amount of RAM to request/use per CPU 



#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
for rep in 1 2 3 4 5 6 7 8 9 10 ; do 
		sbatch --job-name=$jobname \
		--export=JOBNAME=$jobname,SLIMSCRIPT=$slimscript,TREEPROCESS=$treeprocess,KVALUE=$kvalue,BVALUE=$bvalue,REP=$rep,CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir \
		--cpus-per-task=$cpus \
		--mem-per-cpu=$ram_per_cpu \
		--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
		--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
		--time=24:00:00 \
		$executable
done	

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

