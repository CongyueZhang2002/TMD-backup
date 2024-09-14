#! /bin/bash

for i in {0..9}
do
mkdir plots/data_reps/data_rep_00$i
cd -
done

for i in {10..99}
do
mkdir plots/data_reps/data_rep_0$i
cd -

done

for i in {100..200}
do
mkdir plots/data_reps/data_rep_$i
cd -
done
