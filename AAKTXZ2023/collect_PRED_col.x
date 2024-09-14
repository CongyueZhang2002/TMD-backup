#! /bin/bash

# copy results to data_reps directory
mkdir data_reps_PRED_col
for i in {000..106}
do
rm -r data_reps_PRED_col/data_rep_$i/
mkdir data_reps_PRED_col/data_rep_$i/
cp -a replica_$i/plot_data/. data_reps_PRED_col/data_rep_$i/
done
