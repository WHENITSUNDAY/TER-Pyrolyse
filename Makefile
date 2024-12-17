FC=gfortran
FFLAGS = -O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero
EXE=run

all :	$(EXE)

$(EXE) :	constantes.o schema.o main.o
	$(FC) -o $(EXE) constantes.o schema.o main.o

constantes.o :	constantes.f90
	$(FC) $(FFLAGS) -c constantes.f90

schema.o :	schema.f90 constantes.o
	$(FC) $(FFLAGS) -c schema.f90

main.o :	main.f90 constantes.o schema.o
	$(FC) $(FFLAGS) -c main.f90

clean :
	rm -f *.o *.mod $(EXE)