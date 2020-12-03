FC = ifort
EXE = bonCar

objects = main.o

default: $(objects) $(modules)
	ifort -o $(EXE) $(objects) $(modules)
	time ./$(EXE)

module.o: module.f90
	ifort -c module.f90

utility.o: utility.f90
	ifort -c utility.f90

modules = module.o utility.o

main.o: main.f90 $(modules)
	ifort -c main.f90
