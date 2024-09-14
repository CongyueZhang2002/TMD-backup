#!/bin/bash

for i in {100..300}
do
mkdir replica_$i
mkdir replica_$i/weights
mkdir replica_$i/plot_data/
mkdir replica_$i/expdata/
mkdir replica_$i/tools

mkdir replica_$i/plot_data/HERMES/
mkdir replica_$i/plot_data/JLAB2022/
mkdir replica_$i/plot_data/TMD/
mkdir replica_$i/plot_data/PRED/

mkdir replica_$i/numerical 
mkdir replica_$i/TMDs
mkdir replica_$i/TMDs/DSS14_HESSIAN/
mkdir replica_$i/TMDs/DSS17_HESSIAN/ 
mkdir replica_$i/processes
mkdir replica_$i/grid_data

cp fit/Makefile replica_$i/Makefile

cp fit/predict.f replica_$i/predict.f
cp fit/nFF_pre.f replica_$i/nFF_pre.f

cp fit/tools/readdata.f replica_$i/tools/readdata.f
cp fit/tools/data-inc.f replica_$i/tools/data-inc.f

cp fit/best-params.f      replica_$i/best-params.f
cp fit/best-params-rand.f replica_$i/best-params-rand.f

cp fit/submit_predict.sh replica_$i/submit_predict.sh

sed -i "17s#.*#      cd ../replica_$i#" replica_$i/submit_predict.sh

cp -r fit/expdata replica_$i
done

counter=1
j=100

# Read array but skip the first line (header) using awk 'NR>1'
readarray -t my_array < <(awk 'NR>1 {print $13}' params.dat)

for i in "${my_array[@]}"; do
  
  # Get the values starting from second line
  best_values=$(sed -n "$((counter+1))p" params.dat)
  
  read -a arr <<< "$best_values"
  
  # Use the values to create a file with the parameters
  cat << EOF > "replica_${j}/best-params-rand.f"
	gammaBEST = ${arr[0]}
	g3fBEST   = ${arr[1]}
	g3DBEST   = ${arr[2]}
	Nq1BEST   = ${arr[3]}
	gq1BEST   = ${arr[4]}
	dq1BEST   = ${arr[5]}
	Nq2BEST   = ${arr[6]}
	gq2BEST   = ${arr[7]}
	dq2BEST   = ${arr[8]}
	p_10BEST  = ${arr[9]}
	p_11BEST  = ${arr[10]}
        p_12BEST  = ${arr[11]}

	VSTRT(1)  =  gammaBEST
	VSTRT(2)  =  g3fBEST
	VSTRT(3)  =  g3DBEST
	VSTRT(4)  =  Nq1BEST
	VSTRT(5)  =  gq1BEST
	VSTRT(6)  =  dq1BEST
	VSTRT(7)  =  Nq2BEST
	VSTRT(8)  =  gq2BEST
	VSTRT(9)  =  dq2BEST
	VSTRT(10) =  p_10BEST
	VSTRT(11) =  p_11BEST
        VSTRT(12) =  p_12BEST
EOF

  ((counter++))
  ((j++))

done
