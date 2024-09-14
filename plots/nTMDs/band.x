./make_reps.x
./run_replicas.x

baseout=`cat baseline.txt`

myjobs > out.txt
sed -n '2p' out.txt > newout.txt
jobsout=`cat newout.txt`

while [ "$baseout" != "$jobsout" ]
do
    sleep 5

    myjobs > out.txt
    sed -n '2p' out.txt > newout.txt
    jobsout=`cat newout.txt`
enddo

./collect_reps.x

