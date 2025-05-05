FC=gfortran
#FFLAGS = -O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero
FFLAGS = -O3 -fopenmp -march=native

EXE=run

all :	$(EXE)

$(EXE) :	mod_param.o mod_algebre.o mod_fonctions.o mod_schemas.o couplage.o
	$(FC) $(FFLAGS) -o $(EXE) mod_param.o mod_algebre.o mod_fonctions.o mod_schemas.o couplage.o

mod_param.o : mod_param.f90
	$(FC) $(FFLAGS) -c mod_param.f90
	
mod_algebre.o :	mod_algebre.f90 mod_param.o 
	$(FC) $(FFLAGS) -c mod_algebre.f90

mod_fonctions.o : mod_fonctions.f90 mod_param.o
	$(FC) $(FFLAGS) -c mod_fonctions.f90

mod_schemas.o : mod_schemas.f90 mod_param.o mod_algebre.o mod_fonctions.o
	$(FC) $(FFLAGS) -c mod_schemas.f90

couplage.o :	couplage.f90 mod_param.o mod_algebre.o mod_fonctions.o mod_schemas.o
	$(FC) $(FFLAGS) -c couplage.f90

clean :
	rm -f *.o *.mod $(EXE)