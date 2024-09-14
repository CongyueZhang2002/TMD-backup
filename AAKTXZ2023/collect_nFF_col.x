#! /bin/bash

for i in 000 {107..120}
do
cp -a replica_$i/plot_data/TMD/. ../plots/data_reps_nFF_col/data_rep_$i/
done
