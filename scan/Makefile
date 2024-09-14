COMP = gfortran
FFLAGS = -O2	# This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow
SOURCE = fit.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/TMDFF.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS16.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
SRCSCN = scan.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/TMDFF.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS16.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
SRCSC2 = scan2.f ../fit/numerical/gauss.f ../fit/tools/cuts.f90 ../fit/rand.f ../fit/TMDs/TMDFF.f ../fit/TMDs/TMDPDF.f ../fit/TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f ../fit/TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f tools/readdata.f ../fit/numerical/ogata.f ../fit/processes/sidis.f ../fit/processes/dy.f ../fit/TMDs/evolve.f ../fit/TMDs/LIKEn.f  ../fit/TMDs/EPPS16.f ../fit/tools/timer.f ../fit/TMDs/kernel_tmd.f ../fit/numerical/quadpack.f
OBJ=../minuit
LHAPDFLIBS = `lhapdf-config --ldflags`

SOURCE2 = DIS_main.f tools/readdata_DIS.f numerical/gauss.f processes/dis.f processes/sidis.f TMDs/evolve.f TMDs/TMDPDF.f TMDs/TMDFF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f TMDs/kernel_tmd.f numerical/ogata.f TMDs/LIKEn.f TMDs/EPPS16.f
OBJ=../minuit
APFELLIBSMA  = -L/u/home/m/misharya/Projects/apfel/../../lib -lAPFEL
APFELLIBSJT  = -L/u/home/j/jdterry/backup-jdterry/Projects/apfel/../../lib -lAPFEL

SOURCEEIC = EIC.f numerical/gauss.f TMDs/TMDFF.f TMDs/TMDPDF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f numerical/ogata.f processes/sidis.f TMDs/evolve.f TMDs/LIKEn.f  TMDs/EPPS16.f TMDs/kernel_tmd.f numerical/quadpack.f processes/dis.f



executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o $(OBJ)/*.o -o fit.out $<  $(LHAPDFLIBS)

scan:
	$(COMP) $(FFLAGS) -c $(SRCSCN)
	$(COMP) $(LIBPATH) $(LIB) *.o -o scan.out $<  $(LHAPDFLIBS)

scn2:
	$(COMP) $(FFLAGS) -c $(SRCSC2)
	$(COMP) $(LIBPATH) $(LIB) *.o -o scn2.out $<  $(LHAPDFLIBS)

dis:
	$(COMP) $(FFLAGS) -c $(SOURCE2)
	$(COMP) $(LIBPATH) $(LIB) *.o -o dis.out $<  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT)

EIC:
	$(COMP) $(FFLAGS) -c $(SOURCEEIC)
	$(COMP) $(LIBPATH) $(LIB) *.o -o EIC.out $<  $(LHAPDFLIBS) $(APFELLIBSMA) $(APFELLIBSJT)

clean:
	-rm *.out
	-rm *.o
	-rm *~
	-rm *.mod
	-rm fort*
