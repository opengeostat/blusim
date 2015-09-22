
FC=gfortran
#FFLAGS=-Ofast -Wall -mtune=native -msse4  -funroll-loops --param max-unroll-times=4 -DOPTI=1
FFLAGS=-Ofast -Wall -mtune=native -msse4  -funroll-loops --param max-unroll-times=4
OPENMP=-fopenmp

seq: 
	cd gslib90/gslib; make clean; make gslib.a; cd ../..
	$(FC) $(FFLAGS) -c blusim.fpp 
	$(FC) $(FFLAGS) -I. blusim.o -o blusim.exe gslib90/gslib/gslib.a

par: 
	cd gslib90/gslib; make clean; make gslib.a; cd ../..
	$(FC) $(FFLAGS) $(OPENMP) -c blusim.fpp
	$(FC) $(FFLAGS) $(OPENMP) -I. blusim.o -o blusim.exe gslib90/gslib/gslib.a

clean:
	rm *.exe *.o *.mod
