if [ $USER == 'jdterry' ]
then

rt=~/"backup-jdterry/share/LHAPDF/"

else

rt=~/"share/LHAPDF/"

fi

cd $rt
rm -rf nCTEQ15AU
rm -rf nCTEQ15BE
rm -rf nCTEQ15FE
rm -rf nCTEQ15HE
rm -rf nCTEQ15JCC
rm -rf nCTEQ15JFE
rm -rf nCTEQ15JPB
rm -rf nCTEQ15KR
rm -rf nCTEQ15NE
rm -rf nCTEQ15WW
rm -rf nCTEQ15XE
rm -rf nCTEQ15CA

cd -
cp -r nCTEQ15/nCTEQ15AU  "$rt/nCTEQ15AU"
cp -r nCTEQ15/nCTEQ15BE  "$rt/nCTEQ15BE"
cp -r nCTEQ15/nCTEQ15FE  "$rt/nCTEQ15FE"
cp -r nCTEQ15/nCTEQ15HE  "$rt/nCTEQ15HE"
cp -r nCTEQ15/nCTEQ15JCC "$rt/nCTEQ15JCC"
cp -r nCTEQ15/nCTEQ15JFE "$rt/nCTEQ15JFE"
cp -r nCTEQ15/nCTEQ15JPB "$rt/nCTEQ15JPB"
cp -r nCTEQ15/nCTEQ15KR  "$rt/nCTEQ15KR"
cp -r nCTEQ15/nCTEQ15NE  "$rt/nCTEQ15NE"
cp -r nCTEQ15/nCTEQ15WW  "$rt/nCTEQ15WW"
cp -r nCTEQ15/nCTEQ15XE  "$rt/nCTEQ15XE"
cp -r nCTEQ15/nCTEQ15CA  "$rt/nCTEQ15CA"
