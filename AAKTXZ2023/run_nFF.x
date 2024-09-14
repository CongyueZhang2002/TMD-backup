#!/bin/bash

#sleep 1h

for i in 000 {107..120}
do
cd replica_$i
qsub submit_nFF.sh
cd -
done

for i in {001..106} {121..176}
do
rm -r replica_$i
done
