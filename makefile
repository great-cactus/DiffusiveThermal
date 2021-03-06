FC=gfortran
FCFLAGS=-O0 -Wall -fbacktrace -g -fbounds-check
#FCFLAGS=-O2
OBJS=func.o main.o
.SUFFIXES:.f90
driv: $(OBJS)
	$(FC) -o $@ $(FCFLAGS) $(OBJS)

.f90.o:
	$(FC) -o $@ $(FCFLAGS) -c $<
main.o: func.f90
clean:
	rm -f *.o *.mod driv
