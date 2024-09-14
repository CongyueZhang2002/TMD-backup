#! /bin/bash

for i in {100..199}
do
mkdir rep_$i
mkdir rep_$i/weights
mkdir rep_$i/plot_data/
mkdir rep_$i/expdata/
mkdir rep_$i/tools

mkdir rep_$i/plot_data/HERMES/
mkdir rep_$i/plot_data/JLAB2022/
mkdir rep_$i/expdata/JLAB2022/
mkdir rep_$i/numerical 
mkdir rep_$i/TMDs
mkdir rep_$i/TMDs/DSS14_HESSIAN/
mkdir rep_$i/TMDs/DSS17_HESSIAN/ 
mkdir rep_$i/processes
mkdir rep_$i/grid_data

cp fit/Makefile rep_$i/Makefile

cp fit/grid.f rep_$i/grid.f
cp fit/tools/readdata.f rep_$i/tools/readdata.f
cp fit/tools/data-inc.f rep_$i/tools/data-inc.f

cp fit/best-params.f      rep_$i/best-params.f
cp fit/best-params-rand.f rep_$i/best-params-rand.f

cp fit/submit_grid.sh rep_$i/submit_grid.sh

sed -i "17s#.*#      cd ../rep_$i#" rep_$i/submit_grid.sh
sed -i "331s#.*#      j = $i#" rep_$i/grid.f

cp -r fit/expdata rep_$i
done


