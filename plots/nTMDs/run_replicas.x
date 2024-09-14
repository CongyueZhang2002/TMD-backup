for i in {0..9}
do
cd ./replica_00$i/
qsub -cwd -V -N rep_00$i -l h_data=2G,h_rt=50:00:00,highp $PWD/run_fit.bash
cd -
done

for i in {10..99}
do
cd ./replica_0$i/
qsub -cwd -V -N rep_0$i  -l h_data=2G,h_rt=50:00:00,highp $PWD/run_fit.bash
cd -
done

for i in {100..199}
do
cd ./replica_$i/
qsub -cwd -V -N rep_$i   -l h_data=2G,h_rt=50:00:00,highp $PWD/run_fit.bash
cd -
done

cd ./replica_200/
qsub -cwd -V -N rep_200   -l h_data=2G,h_rt=50:00:00,highp -M $USER -m bea  $PWD/run_fit.bash
