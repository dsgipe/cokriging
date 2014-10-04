COMP = g++
FILES = CK.cpp cokriging.o lpkinterface.o vector2array.o array.o
PATH2LAPACK = ../../lapack-3.3.0 
PATH2ARRAY = ../../array/
LAPK = -llapack -lblas
FC = -lgfortran
OUT = CK
INCLUDES = -I $(PATH2ARRAY) -I .
ck: $(FILES)
	$(COMP) $(FILES) -L $(PATH2LAPACK)  $(INCLUDES) $(LAPK) $(FC) -o $(OUT)
cokriging.o : cokriging.cpp cokriging.h
	$(COMP) -g -c $(INCLUDES) cokriging.cpp 
lpkinterface.o : lpkinterface.cpp lpkinterface.h
	$(COMP) -g -c $(INCLUDES) lpkinterface.cpp
vector2array.o : vector2array.cpp vector2array.h
	$(COMP) -g -c vector2array.cpp
array.o : $(PATH2ARRAY)array.cpp $(PATH2ARRAY)array.h
	$(COMP) -g -c $(INCLUDES) $(PATH2ARRAY)array.cpp
clean:
	\rm *.o *~ $(OUT)
