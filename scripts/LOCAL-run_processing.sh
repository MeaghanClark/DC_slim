#!/bin/bash
		
# Last updated 04/28/2022 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

echo "before calling source: $PATH"
source activate slim
echo "after calling source: $PATH"

#load programs we want to use
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module list

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
		python -u ${pythondir}/${treeprocess} ${filename} ${outdir} ${model}_${r}_${rep}_${date} ${mu} ${gen}
	

	done
done	

echo ----------------------------------------------------------------------------------------
echo My executable is $executable		
	

