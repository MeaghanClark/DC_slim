#!/bin/bash

# Last updated 01/13/2023 by MI Clark, script format by R Toczydlowski 

#  run from project directory (where you want output directory to be created)

#define variables:
date=$(date +%m%d%Y)
storagenode=/mnt/home/clarkm89 #path to top level of dir where input/output files live
jobname=load_theta #label for SLURM book-keeping 
run_name=DC_slim #label to use on output files

# define directories
logfilesdir=R_log_theta_${date} #name of directory to create and then write log files to
homedir=$storagenode/$run_name
indir=theta_output_bins_11082022/
outdir=theta_data_obj_${date}

executable=$storagenode/$run_name/scripts/load_theta_into_R.sbatch #script to run 
rscript=load_into.R #R script called by executable


cpus=1 #number of CPUs to request/use per dataset 
ram_per_cpu=75G #amount of RAM to request/use per CPU; 6 G 


#---------------------------------------------------------
#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
#if [ ! -d ./$outdir ]; then mkdir ./$outdir; fi


#submit job to cluster
sbatch --job-name=$jobname \
--export=JOBNAME=$jobname,LOGFILESDIR=$logfilesdir,INDIR=$indir,OUTDIR=$outdir,HOMEDIR=$homedir,RSCRIPT=$rscript \
--cpus-per-task=$cpus \
--mem-per-cpu=$ram_per_cpu \
--output=./$logfilesdir/${jobname}_%A.out \
--error=./$logfilesdir/${jobname}_%A.err \
--time=168:00:00 \
$executable

echo ----------------------------------------------------------------------------------------
echo I am exporting: $jobname, $cpus, $run_name, $storagenode, $outdir, $logfilesdir, $rscript, $datapath, $kvalue, $spatial
echo My executable is $executable               
