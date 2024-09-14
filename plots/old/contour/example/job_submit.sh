#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=50:00:00,h_data=4G,highp
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module use /u/project/CCN/apps/modulefiles

# Load the FSL module
module load fsl

# This is optional
# More info here: https://www.ccn.ucla.edu/wiki/index.php/Hoffman2:FSL
export NO_FSL_JOBS=true

# Your script content goes here...
cd ../../../
source setup.bash
cd plots/contour/example
make clean
make
./example.out
