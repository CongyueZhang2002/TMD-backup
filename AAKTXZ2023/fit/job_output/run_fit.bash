#! /bin/bash
#$ -j y
#$ -pe shared 2

. /u/local/Modules/default/init/modules.sh
module use /u/project/CCN/apps/modulefiles

cd ../../..
source cz.bash
make clean
make 

./fit.out

#cd job_output/
#module load python/2.7.13
#python plot.py
