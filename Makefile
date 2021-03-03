FC = gfortran
FCFLAGS =
LDFLAGS = -llapack
PROJDIR := $(realpath $(CURDIR))
SOURCEDIR := $(PROJDIR)/src

#Objects
sem_forces.o: $(SOURCEDIR)/sem_forces.f90
	$(FC) $(FCFLAGS) -c $(SOURCEDIR)/sem_forces.f90

all: sem_forces.o
	$(FC) $(FCFLAGS) sem_forces.o -o sem_forces $(LDFLAGS)

clean:
	rm -rf *.mod *.o

realclean:
	rm -rf *.mod *.o sem_forces