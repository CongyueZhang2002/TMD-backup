#!/bin/bash

#sleep 1h

for i in {100..300}
do
cd replica_$i
qsub submit_predict.sh
cd -
done
