#! /bin/bash

for i in {0..9}
do
cp fit/predict.f replica_00$i/fit.f
cp fit/best-params.f replica_00$i/best-params.f
done

for i in {10..99}
do
cp fit/predict.f replica_0$i/fit.f
cp fit/best-params.f replica_0$i/best-params.f
done

for i in {100..200}
do
cp fit/predict.f replica_$i/fit.f
cp fit/best-params.f replica_$i/best-params.f
done
