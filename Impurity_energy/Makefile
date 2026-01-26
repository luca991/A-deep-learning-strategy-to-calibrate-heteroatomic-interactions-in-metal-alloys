CC     = g++
CFLAGS = -O2 -march=native -Wall
CFLGS = -lgsl -lgslcblas -lm

all: impure

data.o: data.cpp
	$(CC) $(CFLAGS) -c data.cpp $(CFLGS)

functions.o: functions.cpp
	$(CC) $(CFLAGS) -c functions.cpp $(CFLGS)

impure: impure_energy_calc.cpp data.o functions.o
	$(CC) $(CFLAGS) -o impure impure_energy_calc.cpp data.o functions.o $(CFLGS)


clear:
	rm impure *.o
