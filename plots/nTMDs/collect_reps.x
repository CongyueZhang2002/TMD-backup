#! /bin/bash

# copy results to data_reps directory
for i in {0..9}
do
cp replica_00$i/fort.*   plots/data_reps/data_rep_00$i/
done

for i in {10..99}
do
cp replica_0$i/fort.*  plots/data_reps/data_rep_0$i/
done

for i in {100..200}
do
cp replica_$i/fort.*   plots/data_reps/data_rep_$i/
done
