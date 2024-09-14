#! /bin/bash

make clean
make
./generate.out

mv EPPS16nlo_CT14nlo_He4_0000.dat ./EPPS16nlo_CT14nlo_He4/EPPS16nlo_CT14nlo_He4_0000.dat
mv EPPS16nlo_CT14nlo_Ne20_0000.dat ./EPPS16nlo_CT14nlo_Ne20/EPPS16nlo_CT14nlo_Ne20_0000.dat
mv EPPS16nlo_CT14nlo_Kr84_0000.dat ./EPPS16nlo_CT14nlo_Kr84/EPPS16nlo_CT14nlo_Kr84_0000.dat
mv EPPS16nlo_CT14nlo_Xe131_0000.dat ./EPPS16nlo_CT14nlo_Xe131/EPPS16nlo_CT14nlo_Xe131_0000.dat

cp ./info_files/EPPS16nlo_CT14nlo_He4.info ./EPPS16nlo_CT14nlo_He4/EPPS16nlo_CT14nlo_He4.info
cp ./info_files/EPPS16nlo_CT14nlo_Ne20.info ./EPPS16nlo_CT14nlo_Ne20/EPPS16nlo_CT14nlo_Ne20.info
cp ./info_files/EPPS16nlo_CT14nlo_Kr84.info ./EPPS16nlo_CT14nlo_Kr84/EPPS16nlo_CT14nlo_Kr84.info
cp ./info_files/EPPS16nlo_CT14nlo_Xe131.info ./EPPS16nlo_CT14nlo_Xe131/EPPS16nlo_CT14nlo_Xe131.info

echo "Dont forget to move files to the LHAPDF PDF Set Directory on your system!"
