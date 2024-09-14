#! /bin/bash

for i in {100..300}
do

rm -r data_reps_nFF_pre/data_rep_$i
mkdir data_reps_nFF_pre/data_rep_$i
cp -a replica_$i/plot_data/TMD/.  data_reps_nFF_pre/data_rep_$i/
done

