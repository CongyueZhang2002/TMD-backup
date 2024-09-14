#! /bin/bash

cd ../../LHAPDF-6.2.1/

module load gcc/4.9.3
./configure --prefix=$PWD/..
cd ../bin/
export PATH=$PATH:$PWD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
cd -

module load python/2.7.18
cd ../lib/
export PATH=$PATH:$PWD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
export PYTHONPATH=$PYTHONPATH:$PWD/python3.7.0/site-packages/

cd ../Projects/Nuclear-TMDs/
export PYTHONPATH=$PYTHONPATH:$PWD

cd AAKTXZ2023/QCDNUM/qcdnum-18-00-00
./configure --prefix=/u/project/zkang/misharya/

export PATH=$PATH:/u/project/zkang/misharya/bin 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/project/zkang/misharya/lib
export LD_RUN_PATH=$LD_RUN_PATH:/u/project/zkang/misharya/lib/

cd ../../../