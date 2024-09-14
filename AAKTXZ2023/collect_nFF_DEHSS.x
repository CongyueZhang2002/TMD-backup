#! /bin/bash

for i in 000 {121..176}
do
cp -a replica_$i/plot_data/TMD/. data_reps_nFF_DEHSS/data_rep_$i/
done
