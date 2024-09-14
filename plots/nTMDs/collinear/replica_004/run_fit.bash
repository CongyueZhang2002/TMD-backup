#! /bin/bash

make clean -f Makefile
make -f $PWD/Makefile
./example.out

#cd job_output/
#module load python/2.7.13
#python plot.py
