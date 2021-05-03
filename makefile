FC=gfortran
FCFLAGS=-O2
OBJS=func.o main.o
.SUFFIXES:.f90
drive: $(OBJS)
	$(FC) -o $@ $(FCFLAGS) $(OBJS)

.f90.o:
	$(FC) -c $< -o $@
main.o: func.f90
clean:
	rm -f *.o *.mod drive
