#! /bin/bash

for i in {100..299}
do

mkdir E772_grid/rep_$i/

cp rep_$i/grid_data/E772_C1.dat E772_grid/rep_$i/E772_C1.dat
cp rep_$i/grid_data/E772_C2.dat E772_grid/rep_$i/E772_C2.dat
cp rep_$i/grid_data/E772_C3.dat E772_grid/rep_$i/E772_C3.dat
cp rep_$i/grid_data/E772_C4.dat E772_grid/rep_$i/E772_C4.dat

cp rep_$i/grid_data/E772_Ca1.dat E772_grid/rep_$i/E772_Ca1.dat
cp rep_$i/grid_data/E772_Ca2.dat E772_grid/rep_$i/E772_Ca2.dat
cp rep_$i/grid_data/E772_Ca3.dat E772_grid/rep_$i/E772_Ca3.dat
cp rep_$i/grid_data/E772_Ca4.dat E772_grid/rep_$i/E772_Ca4.dat

cp rep_$i/grid_data/E772_Fe1.dat E772_grid/rep_$i/E772_Fe1.dat
cp rep_$i/grid_data/E772_Fe2.dat E772_grid/rep_$i/E772_Fe2.dat
cp rep_$i/grid_data/E772_Fe3.dat E772_grid/rep_$i/E772_Fe3.dat
cp rep_$i/grid_data/E772_Fe4.dat E772_grid/rep_$i/E772_Fe4.dat

cp rep_$i/grid_data/E772_W1.dat E772_grid/rep_$i/E772_W1.dat
cp rep_$i/grid_data/E772_W2.dat E772_grid/rep_$i/E772_W2.dat
cp rep_$i/grid_data/E772_W3.dat E772_grid/rep_$i/E772_W3.dat
cp rep_$i/grid_data/E772_W4.dat E772_grid/rep_$i/E772_W4.dat 

done
