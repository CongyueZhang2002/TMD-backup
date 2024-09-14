#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 1
#$ -l h_rt=1:00:00,h_data=4G,highp
#
#

. /u/local/Modules/default/init/modules.sh
module use /u/project/CCN/apps/modulefiles


# Your script content goes here...
if [ "$USER" == "maxzhang" ]; then
    cd ../../
    source cz.bash
    cd ../fit
    make clean
    make c
    ./fitc.out
fi
#if [ "$USER" == "" ]; then
    #cd ../../
    #chmod +x setup.bash
    #source setup.bash
    #chmod +x Makefile
    #make clean
    #make
    #chmod +x fit.out
    #./fit.out
#fi
