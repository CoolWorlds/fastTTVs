NESTLIBDIR = ..
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
FC=mpif90
OBJFILES = params.o inversebeta.o jasmine.o mandel.o plan.o transit.o radial.o like.o nestwrap.o main.o

all: LRVplan 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
LRVplan: $(OBJFILES)
	$(FC) -o ../LRVplan $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../LRVplan
