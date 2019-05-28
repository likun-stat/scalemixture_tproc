
#!/bin/bash

#PBS -A bas59_a_g_bc_default
#PBS -V
#PBS -l nodes=1:ppn=20
#PBS -l walltime=48:00:00
#PBS -N tproc_DA_cont
#PBS -j oe
#PBS -l pmem=3gb

#-#-# Echo
echo "#-#-#Job started on `hostname` at `date` "
echo This job runs on the following processors:
echo `cat $PBS_NODEFILE`


cd $PBS_O_WORKDIR
echo " "
echo " "
echo "Job started on hostname at date"


module purge
module load r/3.4

Rscript imputation_animation_02_cont.R

echo "#-#-#Job Ended at `date`"
