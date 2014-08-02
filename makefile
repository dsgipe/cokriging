COMP = g++
FILES = CK.cpp cokriging.o lpkinterface.o vector2array.o
PATH2LAPACK = ../lapack-3.3.0 
LAPK = -llapack -lblas
FC = -lgfortran
OUT = CK

ck: $(FILES)
	$(COMP) $(FILES) -L $(PATH2LAPACK) $(LAPK) $(FC) -o $(OUT)
cokriging.o : cokriging.cpp cokriging.h
	$(COMP) -g -c cokriging.cpp
lpkinterface.o : lpkinterface.cpp lpkinterface.h
	$(COMP) -g -c lpkinterface.cpp
vector2array.o : vector2array.cpp vector2array.h
	$(COMP) -g -c vector2array.cpp
clean:
	\rm *.o *~ $(OUT)
