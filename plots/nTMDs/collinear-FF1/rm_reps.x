#! /bin/bash

#cd ../
#chmod -R 777 AAKTX2021
#cd AAKTX2021

for i in {0..9}
do
rm -rf replica_00$i
done

for i in {10..99}
do
rm -rf replica_0$i
done

for i in {100..200}
do
rm -rf replica_$i
done
