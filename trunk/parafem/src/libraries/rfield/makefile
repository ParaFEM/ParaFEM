
# Source code version = 5.39

FC      = ftn
FFLAGS  = -V -O3 -s real64

FILES =	sim3de.f

OBJS  =	sim3de.o

.SUFFIXES: .f .o

.f.o:
	$(FC) -c $(FFLAGS) *.f

lib: $(OBJS)
	ar -r rfield.a $(OBJS)
	cp rfield.a ../../../bin
	
clean:
	rm *.o
