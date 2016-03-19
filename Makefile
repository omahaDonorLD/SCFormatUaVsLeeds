################################################################
#                  VARIABLES DEFINITIONS                       #
################################################################

# Compilator
CPP = g++

# Options of compilation
CFLAGS = -Wall -Wextra -pedantic -ansi

# Options of compilation for performances' optimisation
OPTIMIZATION = -O3

# Folder who contained files to include
INCLUDE = -I .

# Finals executables
AIMS = NSGA2

################################################################
#                  COMMANDS OF COMPILATION                     #
################################################################

all: $(AIMS)

NSGA2: individu.o NSGA-2.o individuFactory.o output.o
	$(CPP) -o NSGA2 individu.o NSGA-2.o individuFactory.o output.o $(CFLAGS) $(OPTIMIZATION)

individuFactory.o: individu.hpp individuFactory.cpp
	$(CPP) -o individuFactory.o -c individuFactory.cpp $(CFLAGS) $(OPTIMIZATION)

individu.o: individu.cpp
	$(CPP) -o individu.o -c individu.cpp $(CFLAGS) $(OPTIMIZATION)
	
output.o: individu.hpp output.cpp
	$(CPP) -o output.o -c output.cpp $(CFLAGS) $(OPTIMIZATION)

NSGA-2.o: NSGA-2.cpp individu.hpp individuFactory.hpp output.hpp
	$(CPP) -o NSGA-2.o -c NSGA-2.cpp $(CFLAGS) $(OPTIMIZATION)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(AIMS)
	
blank:
	tput clear
