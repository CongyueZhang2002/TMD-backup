This directory analyzes the Drell Yan dilepton production data.

Data is read in readdata.f
fit.f is the main program.

To run the fit, use
make
./fit.out

To run a prediction use
make predict 
./predict.out

To make predictions for EIC use
make EIC
./EIC.out


note that the scheme can be changed by editing the fit.f code

LO -> nloops = 1  (kind of a misnomer)
NLO-> nloops = 2

LL  -> nll = 1
NLL -> nll = 2
NNLL-> nll = 3

ogata -> fft = 1
gauss -> fft = 0

hoppet -> hop = 1 (Very slow and only moderately better for describing the data)
dglap  -> hop = 0

Plotting software is done in python in the plots directory.
For use in hoffman use python/2.7.13, may move to python 3.7.0.
