#!/bin/bash

#sleep 1h

for i in {101..300}
do
cd replica_$i
qsub submit2.sh
cd -
done

