#! /bin/bash

for i in {0..9}
do
find ./replica_00$i -type f -name "*rep_00$i*" -exec rm {} \;
done

for i in {10..99}
do
find ./replica_0$i -type f -name "*rep_0$i*" -exec rm {} \;
done

for i in {100..200}
do
find ./replica_$i -type f -name "*rep_$i*" -exec rm {} \;
done

for i in {0..9}
do
sed -i "24s/.*/    cd replica_00$i/" replica_00$i/job_submit.sh
sed -i "5s/.*/#$ -pe shared 2/" replica_00$i/job_submit.sh
rm replica_00$i/fit.f
rm replica_00$i/best-params-rand.f
rm replica_00$i/result.txt
cp fit/fit.f replica_00$i
cp fit/best-params-rand.f replica_00$i
done

for i in {10..99}
do
sed -i "24s/.*/    cd replica_0$i/" replica_0$i/job_submit.sh
sed -i "5s/.*/#$ -pe shared 2/" replica_0$i/job_submit.sh
rm replica_0$i/fit.f
rm replica_0$i/result.txt
rm replica_0$i/best-params-rand.f
cp fit/fit.f replica_0$i
cp fit/best-params-rand.f replica_0$i
done

for i in {100..200}
do
sed -i "24s/.*/    cd replica_$i/" replica_$i/job_submit.sh
sed -i "5s/.*/#$ -pe shared 2/" replica_$i/job_submit.sh
rm replica_$i/fit.f
rm replica_$i/result.txt
rm replica_$i/best-params-rand.f
cp fit/fit.f replica_$i
cp fit/best-params-rand.f replica_$i
done
