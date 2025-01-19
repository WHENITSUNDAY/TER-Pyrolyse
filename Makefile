FC=gfortran
#FFLAGS = -O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero
FFLAGS = -O3 -march=native
EXE=run

all :	$(EXE)

$(EXE) :	constantes.o schema.o chimie.o
	$(FC) -o $(EXE) constantes.o schema.o chimie.o

constantes.o :	constantes.f90
	$(FC) $(FFLAGS) -c constantes.f90

schema.o :	schema.f90 constantes.o
	$(FC) $(FFLAGS) -c schema.f90

chimie.o :	chimie.f90 constantes.o schema.o
	$(FC) $(FFLAGS) -c chimie.f90

clean :
	rm -f *.o *.mod $(EXE)