#! /bin/bash

# copy results to data_reps directory
for i in {000..106}
do
cp -a replica_$i/plot_data/. ../plots/data_reps_col/data_rep_$i/
done


