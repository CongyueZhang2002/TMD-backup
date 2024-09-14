COMP = gfortran
FFlAGS = -O
SOURCE = example.f ../example/TMDPDF.f ../example/gauss.f ../example/ogata.f ../example/kernel_tmd.f ../example/quadpack.f ../example/common.f readdata.f ../example/sivers.f ../example/TMDFF.f ../example/DSS14_HESSIAN/fDSS_HESSIAN14.f ../example/DSS17_HESSIAN/fDSS_HESSIAN17.f ../example/fDSS.f ../example/FF_NLO.f ../example/EPPS16.f ../example/PDF_NLO.f
LHAPDFLIBS = `lhapdf-config --ldflags`

example:
	$(COMP) $(FFLAGS) -c $(SOURCE)

	$(COMP) *.o -o example.out $<  $(LHAPDFLIBS)

clean:
	rm *.o
	rm *.out
