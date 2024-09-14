#! /bin/bash

for i in {100..199}
do

mkdir pwr_grid/rep_$i
mkdir pwr_grid/rep_$i/grid_data

for j in {1..4}
do
cp rep_$i/grid_data/E772_C$j.dat  pwr_grid/rep_$i/grid_data
cp rep_$i/grid_data/E772_Ca$j.dat  pwr_grid/rep_$i/grid_data
cp rep_$i/grid_data/E772_Fe$j.dat  pwr_grid/rep_$i/grid_data
cp rep_$i/grid_data/E772_W$j.dat  pwr_grid/rep_$i/grid_data
done

for j in {1..8}
do
cp rep_$i/grid_data/CMS5_$j.dat  pwr_grid/rep_$i/grid_data
done

for j in {1..7}
do
cp rep_$i/grid_data/ATLAS1_$j.dat  pwr_grid/rep_$i/grid_data
done

for j in {1..4}
do
cp rep_$i/grid_data/RHIC_$j.dat  pwr_grid/rep_$i/grid_data
done

for j in #{1..29}
do
cp rep_$i/grid_data/E866_Fe$j.dat  pwr_grid/rep_$i/grid_data
cp rep_$i/grid_data/E866_W$j.dat  pwr_grid/rep_$i/grid_data
done

done
