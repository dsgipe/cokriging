COMP = g++
FILES = CK.cpp cokriging.h cokriging.cpp lpkinterface.h vector2array.h
PATH2LAPACK = ../lapack-3.3.0 
LAPK = -llapack -lblas
FC = -lgfortran
OUT = CK

ck: $(FILES)
	$(COMP) $(FILES) -L $(PATH2LAPACK) $(LAPK) $(FC) -o $(OUT)

clean:
	\rm *.o *~ $(OUT)
