#!/bin/bash

#sleep 1h

for i in {000..106}
do
cd replica_$i
qsub submit_c.sh
cd -
done

for i in {107..176}
do
rm -r replica_$i
done
