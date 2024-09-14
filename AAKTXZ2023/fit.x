#! /bin/bash

for i in {303..306}
do
mkdir replica_$i
mkdir replica_$i/weights
mkdir replica_$i/plot_data/
mkdir replica_$i/expdata/
mkdir replica_$i/tools

mkdir replica_$i/expdata/DRELLYAN/
mkdir replica_$i/expdata/HERMES_DIS/
mkdir replica_$i/expdata/DRELLYAN/E288/
mkdir replica_$i/expdata/DRELLYAN/E605/
mkdir replica_$i/expdata/DRELLYAN/E772/
mkdir replica_$i/expdata/DRELLYAN/E866/
mkdir replica_$i/expdata/DRELLYAN/CMS8/
mkdir replica_$i/expdata/DRELLYAN/CMS5/
mkdir replica_$i/expdata/DRELLYAN/ATLAS3/
mkdir replica_$i/expdata/DRELLYAN/ATLAS5/
mkdir replica_$i/expdata/DRELLYAN/E906/
mkdir replica_$i/expdata/DRELLYAN/RHIC/
mkdir replica_$i/expdata/DRELLYAN/RHIC2/
mkdir replica_$i/expdata/SIDIS/
mkdir replica_$i/expdata/SIDIS/HERMES_DIS/
mkdir replica_$i/expdata/SIDIS/JLAB_DIS/
mkdir replica_$i/expdata/SIDIS/JLAB12_DIS/
mkdir replica_$i/expdata/JLAB2022/

mkdir replica_$i/plot_data/HERMES/
mkdir replica_$i/plot_data/JLAB2022/
mkdir replica_$i/plot_data/TMD/

mkdir replica_$i/numerical 
mkdir replica_$i/TMDs
mkdir replica_$i/TMDs/DSS14_HESSIAN/
mkdir replica_$i/TMDs/DSS17_HESSIAN/ 
mkdir replica_$i/processes
mkdir replica_$i/grid_data

cp fit/expdata/DRELLYAN/RHIC2/RHIC_pp.dat replica_$i/expdata/DRELLYAN/RHIC2/RHIC_pp.dat
cp fit/expdata/JLAB2022/pi+C.dat replica_$i/expdata/JLAB2022
cp fit/expdata/JLAB2022/pi+Fe.dat replica_$i/expdata/JLAB2022
cp fit/expdata/JLAB2022/pi+Pb.dat replica_$i/expdata/JLAB2022
cp fit/expdata/JLAB2022/pi-C.dat replica_$i/expdata/JLAB2022
cp fit/expdata/JLAB2022/pi-Fe.dat replica_$i/expdata/JLAB2022
cp fit/expdata/JLAB2022/pi-Pb.dat replica_$i/expdata/JLAB2022

cp fit/Makefile replica_$i/Makefile

cp fit/fit.f replica_$i/fit.f
cp fit/fit1.f replica_$i/fit1.f
cp fit/fit2.f replica_$i/fit2.f
#cp fit/fit3.f replica_$i/fit3.f
cp fit/predict.f replica_$i/predict.f
cp fit/fitc.f replica_$i/fitc.f

cp fit/tools/readdata.f replica_$i/tools/readdata.f
cp fit/tools/data-inc.f replica_$i/tools/data-inc.f

cp fit/best-params.f      replica_$i/best-params.f
cp fit/best-params-rand.f replica_$i/best-params-rand.f

cp fit/submit.sh replica_$i/submit.sh
cp fit/submit1.sh replica_$i/submit1.sh
cp fit/submit2.sh replica_$i/submit2.sh
cp fit/submit3.sh replica_$i/submit3.sh

sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit.sh
sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit1.sh
sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit2.sh
sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit3.sh

#m=$(echo "0 + 0.04 * ($i - 100)" | bc)
#sed -i "47s#.*#      NG2BEST = $m#" replica_$i/best-params-rand.f
#sed -i "83s#.*#      NG1BEST = $m#" replica_$i/best-params-rand.f

cp -r fit/expdata replica_$i
done


