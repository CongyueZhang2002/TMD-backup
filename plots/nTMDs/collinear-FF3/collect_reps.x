#! /bin/bash

# copy results to data_reps directory
for i in {0..9}
do
cp replica_00$i/fort.*   plots/data_reps/data_rep_00$i/
done

for i in {10..56}
do
cp replica_0$i/fort.*  plots/data_reps/data_rep_0$i/
done
