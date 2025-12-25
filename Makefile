CC     = g++
CFLAGS = -O2 -march=native
CFLGS = -lgsl -lgslcblas -lm
CFPAR = -fopenmp

all: Unit_cell

data.o: data.cpp
	$(CC) $(CFLAGS) -c data.cpp
	
functions.o: functions.cpp
	$(CC) $(CFLAGS) -c functions.cpp $(CFLGS)

Unit_cell: Unit_cell.cpp data.o functions.o
	$(CC) $(CFLAGS) -o Unit_cell Unit_cell.cpp data.o functions.o $(CFLGS)

clear:
	rm Unit_cell Unit_cell_short Unit_cell_notintcompo *.o
