#! /bin/bash

for i in {000..176} 
do
mkdir replica_$i
mkdir replica_$i/job_output/
mkdir replica_$i/plot_data/
mkdir replica_$i/expdata/
mkdir replica_$i/tools

mkdir replica_$i/plot_data/HERMES/
mkdir replica_$i/expdata/JLAB2022/
mkdir replica_$i/plot_data/JLAB2022/
mkdir replica_$i/plot_data/TMD/
mkdir replica_$i/plot_data/PRED/

cp fit/tools/readdatac.f replica_$i/tools/readdatac.f
cp fit/tools/data-inc.f replica_$i/tools/data-inc.f
cp fit/Makefile replica_$i/Makefile
cp fit/submit_c.sh      replica_$i/submit_c.sh
cp fit/submit_c2.sh      replica_$i/submit_c2.sh
cp fit/fitc.f replica_$i/fitc.f
cp fit/nFF_col.f replica_$i/nFF_col.f
cp fit/nFF_DEHSS.f replica_$i/nFF_DEHSS.f
cp fit/best-params.f      replica_$i/best-params.f
cp fit/best-params-rand.f replica_$i/best-params-rand.f
cp fit/submit_nFF.sh      replica_$i/submit_nFF.sh

sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit_c.sh
sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit_c2.sh
sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit_nFF.sh

cp -r fit/expdata replica_$i
mkdir replica_$i/weights
done

python cp_results.py
