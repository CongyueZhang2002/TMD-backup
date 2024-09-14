#! /bin/bash

for i in {0..9}
do
rm -rf replica_00$i
done
for i in {10..99}
do
rm -rf replica_0$i
done
for i in {100..400}
do
rm -rf replica_$i
rm -rf rep_$i
done
