COMP = gfortran -g
FFLAGS = -O2    # This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow
SOURCE = trial_QCDNUM.f LIKEn_test.f new_FF_QCDNUM.f DSS14_HESSIAN/fDSS_HESSIAN14.f DSS17_HESSIAN/fDSS_HESSIAN17.f evolve.f 
LHAPDFLIBS = `lhapdf-config --ldflags`
QCDNUMLIBS = `qcdnum-config --ldflags`
APFELLIBSMA  = -L/u/home/m/misharya/Projects/apfel/../../lib -lAPFEL
APFELLIBSJT  = -L/u/home/j/jdterry/backup-jdterry/Projects/apfel/../../lib -lAPFEL
APFELLIBSCZ  = -L/u/home/m/maxzhang/project-zkang/lib -lAPFEL

executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o  -o fit.out $<  $(LHAPDFLIBS) $(QCDNUMLIBS) $(APFELLIBSMA) $(APFELLIBSJT) $(APFELLIBSCZ)

clean:
	-rm *.out
	-rm *.o
	-rm *~
	-rm *.mod
	-rm fort*
	-rm weights/*

