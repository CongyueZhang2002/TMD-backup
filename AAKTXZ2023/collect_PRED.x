#! /bin/bash
mkdir data_reps_PRED
for i in {100..300}
do

rm -r data_reps_PRED/data_rep_$i
mkdir data_reps_PRED/data_rep_$i
cp -a replica_$i/plot_data/.  data_reps_PRED/data_rep_$i/
done
