FC=gfortran
FFLAGS = -O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero
#FFLAGS = -O3 -march=native
EXE=run

all :	$(EXE)

$(EXE) :	mod_constantes.o mod_algebre.o mod_schema.o couplage.o
	$(FC) -o $(EXE) mod_constantes.o mod_algebre.o mod_schema.o couplage.o

mod_constantes.o :	mod_constantes.f90
	$(FC) $(FFLAGS) -c mod_constantes.f90

mod_algebre.o :	mod_algebre.f90 mod_constantes.o
	$(FC) $(FFLAGS) -c mod_algebre.f90

mod_schema.o :	mod_schema.f90 mod_algebre.o mod_constantes.o
	$(FC) $(FFLAGS) -c mod_schema.f90

couplage.o :	couplage.f90  mod_constantes.o mod_algebre.o mod_schema.o
	$(FC) $(FFLAGS) -c couplage.f90

clean :
	rm -f *.o *.mod $(EXE)