#! /bin/bash

if [ $USER == 'jdterry' ]
then

rt=~/"backup-jdterry/share/LHAPDF/"

else

rt=~/"share/LHAPDF/"

fi

cd $rt
rm -rf LIKEnVC
rm -rf LIKEnHE
rm -rf LIKEnNE
rm -rf LIKEnKR
rm -rf LIKEnXE
rm -rf LIKEnCC
rm -rf LIKEnFE
rm -rf LIKEnPB

cd -
cp -r LIKEn/LIKEnVC  "$rt/LIKEnVC"
cp -r LIKEn/LIKEnHE  "$rt/LIKEnHE"
cp -r LIKEn/LIKEnNE  "$rt/LIKEnNE"
cp -r LIKEn/LIKEnKR  "$rt/LIKEnKR"
cp -r LIKEn/LIKEnXE  "$rt/LIKEnXE"
cp -r LIKEn/LIKEnCC    "$rt/LIKEnCC"
cp -r LIKEn/LIKEnFE    "$rt/LIKEnFE"
cp -r LIKEn/LIKEnAU    "$rt/LIKEnAU"
cp -r LIKEn/LIKEnPB    "$rt/LIKEnPB"
