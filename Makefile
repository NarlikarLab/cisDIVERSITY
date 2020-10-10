#CC=gcc -Wall -g -pg -fopenmp -std=c99 -std=gnu99
#CC=gcc -fopenmp -std=c99 -std=gnu99 

CC=gcc -gdwarf-2  -fopenmp -std=c99 -std=gnu99 

LFLAGS=-lm

all: learnDiverseModules 

learnDiverseModules: main.o model.o initializations.o littleHelpers.o  probabilityFunctions.o printing.o samplingFunctions.o
	$(CC) -o learnDiverseModules main.o model.o initializations.o littleHelpers.o probabilityFunctions.o printing.o samplingFunctions.o $(LFLAGS)


clean:
	rm -f *.o
