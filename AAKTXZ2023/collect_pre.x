#! /bin/bash

for i in {100..300}
do

rm -r data_reps/data_rep_$i
mkdir data_reps/data_rep_$i
cp -a replica_$i/plot_data/.  data_reps/data_rep_$i/
done

