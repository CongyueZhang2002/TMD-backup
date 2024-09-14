COMP = gfortran
FFLAGS = -O2	# This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow
SOURCE = example.f LIKEn_test.f DSS14_HESSIAN/fDSS_HESSIAN14.f DSS17_HESSIAN/fDSS_HESSIAN17.f evolve.f
QCDNUMLIBS = `qcdnum-config --ldflags`

executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o -o example.out $<  $(QCDNUMLIBS) 

clean:
	-rm *.out
	-rm *.o
	-rm *~
	-rm *.mod
	-rm fort*
	-rm weights/*
