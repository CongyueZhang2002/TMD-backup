#! /bin/bash

# copy results to data_reps directory
for i in {0..9}
do
cp -a replica_00$i/plot_data/. ../plots/data_reps_col_TMD/data_rep_00$i/
done

for i in {10..96}
do
cp -a replica_0$i/plot_data/. ../plots/data_reps_col_TMD/data_rep_0$i/
done

for i in 
do
cp -a replica_$i/plot_data/. ../plots/data_reps_col_TMD/data_rep_$i/
done
