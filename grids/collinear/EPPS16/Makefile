COMP = gfortran
FFLAGS = -O2	# This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow
SOURCE = generate.f  TMDs/TMDPDF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f numerical/ogata.f TMDs/evolve.f TMDs/EPPS16.f TMDs/kernel_tmd.f numerical/quadpack.f
SOURCE2 = test.f     TMDs/TMDPDF.f TMDs/DSS14_HESSIAN/fDSS_HESSIAN14.f TMDs/DSS17_HESSIAN/fDSS_HESSIAN17.f numerical/ogata.f TMDs/evolve.f TMDs/EPPS16.f TMDs/kernel_tmd.f numerical/quadpack.f
LHAPDFLIBS = `lhapdf-config --ldflags`

executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o -o generate.out $<  $(LHAPDFLIBS)

clean:
	-rm *.out
#	-rm *mod
	-rm *.o
#	-rm fort.*
	-rm *~

test:
	$(COMP) $(FFLAGS) -c $(SOURCE2)
	$(COMP) $(LIBPATH) $(LIB) *.o -o test.out $<  $(LHAPDFLIBS)
