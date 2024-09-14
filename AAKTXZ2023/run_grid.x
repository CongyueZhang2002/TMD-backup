#!/bin/bash

for i in {100..199}
do
cd rep_$i
qsub submit_grid.sh
cd -
done

