#! /bin/bash

if [ $USER == 'jdterry' ]
then

rt=~/"backup-jdterry/share/LHAPDF/"

else

rt=~/"share/LHAPDF/"

fi

cd $rt
rm -rf EPPS16nlo_CT14nlo_He4
rm -rf EPPS16nlo_CT14nlo_Ne20
rm -rf EPPS16nlo_CT14nlo_Kr84
rm -rf EPPS16nlo_CT14nlo_Xe131

cd -
cp -r EPPS16/EPPS16nlo_CT14nlo_He4  "$rt/EPPS16nlo_CT14nlo_He4"
cp -r EPPS16/EPPS16nlo_CT14nlo_Ne20 "$rt/EPPS16nlo_CT14nlo_Ne20"
cp -r EPPS16/EPPS16nlo_CT14nlo_Kr84 "$rt/EPPS16nlo_CT14nlo_Kr84"
cp -r EPPS16/EPPS16nlo_CT14nlo_Xe131 "$rt/EPPS16nlo_CT14nlo_Xe131"
