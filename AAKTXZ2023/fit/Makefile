COMP = gfortran -g
FFLAGS = -O2	# This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow

SOURCE = fit.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
SOURCE1 = fit1.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM1.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis1.f ../fit/processes/dy1.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
SOURCE2 = fit2.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM2.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis2.f ../fit/processes/dy2.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f 
SOURCE3 = fit3.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM3.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis3.f ../fit/processes/dy3.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f 
 
SOURCEPREDICT = predict.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT.f ../fit/TMDs/TMDFF_kT.f
SOURCEPREDICT1 = predict.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM1.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis1.f ../fit/processes/dy1.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT1.f ../fit/TMDs/TMDFF_kT1.f
SOURCEPREDICT2 = predict.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM2.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis2.f ../fit/processes/dy2.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT2.f ../fit/TMDs/TMDFF_kT2.f
SOURCEPREDICT3 = predict.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM3.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis3.f ../fit/processes/dy3.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT3.f ../fit/TMDs/TMDFF_kT3.f

SOURCEGRID = grid.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f 
SOURCEc = fitc.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdatac.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT.f ../fit/TMDs/TMDFF_kT.f
SOURCEc2 = fitc.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM2.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdatac.f ../fit/numerical/ogata.f ../fit/processes/sidis2.f ../fit/processes/dy2.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT2.f ../fit/TMDs/TMDFF_kT2.f
SOURCEPREDICT_TMD = predict_TMD.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f ../fit/TMDs/TMDPDF_kT.f ../fit/TMDs/TMDFF_kT.f
 
SOURCEPREDICTAPFEL = predict_APFEL.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
SRCSCN = scan.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/TMDFF.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
SRCSC2 = scan2.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/TMDFF.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
OBJ=../minuit
LHAPDFLIBS = `lhapdf-config --ldflags`
QCDNUMLIBS = `qcdnum-config --ldflags`

OBJ=../minuit
APFELLIBSMA  = -L/u/home/m/misharya/Projects/apfel/../../lib -lAPFEL
APFELLIBSJT  = -L/u/home/j/jdterry/backup-jdterry/Projects/apfel/../../lib -lAPFEL
APFELLIBSCZ  = -L/u/home/m/maxzhang/project-zkang/lib -lAPFEL

SOURCEEIC = EIC.f numerical/gauss.f TMDs/TMDFF.f TMDs/TMDPDF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f numerical/ogata.f processes/sidis.f TMDs/evolve.f TMDs/LIKEn.f  TMDs/EPPS21.f TMDs/kernel_tmd.f numerical/quadpack.f processes/dis.f
SOURCEFFF = checkFF.f numerical/gauss.f TMDs/TMDPDF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f numerical/ogata.f processes/sidis.f TMDs/evolve.f TMDs/LIKEn.f  TMDs/EPPS21.f TMDs/kernel_tmd.f numerical/quadpack.f processes/dis.f TMDs/new_FF.f
SOURCEDIS = DIS_main.f tools/readdata_DIS.f numerical/gauss.f processes/dis.f processes/sidis.f TMDs/new_FF.f TMDs/evolve.f TMDs/TMDPDF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f TMDs/kernel_tmd.f numerical/ogata.f TMDs/LIKEn.f TMDs/EPPS21.f

Sourcenff_DEHSS = nFF_DEHSS.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f tools/readdatac.f ../fit/numerical/ogata.f  ../fit/TMDs/evolve.f  ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f  
Sourcenff_col = nFF_col.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/TMDFF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdatac.f ../fit/numerical/ogata.f  ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f  
Sourcenff_pre = nFF_pre.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/new_FF_QCDNUM.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS21.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f

executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fit.out $<  $(QCDNUMLIBS)  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
1:
	$(COMP) $(FFLAGS) -c $(SOURCE1)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fit1.out $<  $(QCDNUMLIBS)  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
2:
	$(COMP) $(FFLAGS) -c $(SOURCE2)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fit2.out $<  $(QCDNUMLIBS)  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
3:
	$(COMP) $(FFLAGS) -c $(SOURCE3)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fit3.out $<  $(QCDNUMLIBS)  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

predict:
	$(COMP) $(FFLAGS) -c $(SOURCEPREDICT)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o predict.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
predict1:
	$(COMP) $(FFLAGS) -c $(SOURCEPREDICT1)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o predict.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
predict2:
	$(COMP) $(FFLAGS) -c $(SOURCEPREDICT2)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o predict.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
predict3:
	$(COMP) $(FFLAGS) -c $(SOURCEPREDICT3)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o predict.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

grid:
	$(COMP) $(FFLAGS) -c $(SOURCEGRID)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o grid.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
c:
	$(COMP) $(FFLAGS) -c $(SOURCEc)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fitc.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
c2:
	$(COMP) $(FFLAGS) -c $(SOURCEc2)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fitc.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

predict_TMD:
	$(COMP) $(FFLAGS) -c $(SOURCEPREDICT_TMD)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o predict_TMD.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)



apfel:
	$(COMP) $(FFLAGS) -c $(SOURCEPREDICTAPFEL)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o apfel.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

scan:
	$(COMP) $(FFLAGS) -c $(SRCSCN)
	$(COMP) $(LIBPATH) $(LIB) *.o -o scan.out $<  $(LHAPDFLIBS)

scn2:
	$(COMP) $(FFLAGS) -c $(SRCSC2)
	$(COMP) $(LIBPATH) $(LIB) *.o -o scn2.out $<  $(LHAPDFLIBS)

dis:
	$(COMP) $(FFLAGS) -c $(SOURCEDIS)
	$(COMP) $(LIBPATH) $(LIB) *.o -o dis.out $<  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

EIC:
	$(COMP) $(FFLAGS) -c $(SOURCEEIC)
	$(COMP) $(LIBPATH) $(LIB) *.o -o EIC.out $<  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
checkFF:
	$(COMP) $(FFLAGS) -c $(SOURCEFFF)
	$(COMP) $(LIBPATH) $(LIB) *.o -o checkFF.out $<  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

nFF_col:
	$(COMP) $(FFLAGS) -c $(Sourcenff_col)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o nFF_col.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
nFF_pre:
	$(COMP) $(FFLAGS) -c $(Sourcenff_pre)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o nFF_pre.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)
nFF_DEHSS:
	$(COMP) $(FFLAGS) -c $(Sourcenff_DEHSS)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o nFF_DEHSS.out $< $(QCDNUMLIBS) $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)


clean:
	-rm *.out
	-rm *.o
	-rm *~
	-rm *.mod
	-rm fort*
	-rm weights/*
