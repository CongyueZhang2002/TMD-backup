#! /bin/bash

for i in {000..200}
do
mkdir replica_$i
mkdir replica_$i/expdata/DRELLYAN/
mkdir replica_$i/expdata/DRELLYAN/E906/
#make directories for your experimental datasets

cp fit/Makefile replica_$i/Makefile
#Makefile gives instruction on which files should be compiled together

cp fit/fit.f replica_$i/fit.f
#file for performing fit

cp fit/best-params-rand.f replica_$i/best-params-rand.f 
#set of initial parameters

cp fit/submit.sh replica_$i/submit.sh \
#submit.sh for submittign jobs to Hoffman2

mkdir replica_$i/numerical 
#Other directories and files required for you fit to run

done

python copy_data.py 
#Add gaussian noise to replica data
