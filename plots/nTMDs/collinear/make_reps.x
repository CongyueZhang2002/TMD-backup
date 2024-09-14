#! /bin/bash

for i in {0..9}
do
mkdir replica_00$i
cp example/readdata0.f replica_00$i/readdata.f
cp example/example.f replica_00$i/example.f
cp example/Makefile replica_00$i/Makefile
cp example/run_fit.bash replica_00$i/run_fit.bash
done

for i in {10..96}
do
mkdir replica_0$i
cp example/readdata0.f replica_0$i/readdata.f
cp example/example.f replica_0$i/example.f
cp example/Makefile replica_0$i/Makefile
cp example/run_fit.bash replica_0$i/run_fit.bash
done



#python copy_data.py
#module load python/2.7.13
python cp_results.py
# python copy_readdatas.py
