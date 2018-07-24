################################################################
#                  VARIABLES DEFINITIONS                       #
################################################################

# Compilator
CPP = g++

# Directories
SRCS = srcs
HEADS = heads
OBJS = objs

# include standard math library
LDFLAGS += -lm

# Options of compilation
CFLAGS = -Wall -Wextra -pedantic -ansi

# Options of compilation for performances' optimisation
OPTI = -O3
DBG=-g

# Folder who contained files to include
INCLUDE = -I .

# Finals executables
AIMS = NSGA2

# Igraph libraries
# Personal computer
IGRAPH_INCLUDE=/home/ubuntu/repos_git/opt/igraph-0.7.1/include/
IGRAPH_LIB=/home/ubuntu/repos_git/opt/igraph-0.7.1/lib/
# Auxiliary computer
2ND_IGRAPH_INCLUDE=/home/ubuntu/repos_git/opt/igraph-0.7.1/include/
2ND_IGRAPH_INCLUDE=/home/ubuntu/repos_git/opt/igraph-0.7.1/lib/

################################################################
#                  COMMANDS OF COMPILATION                     #
################################################################

all: $(AIMS)

NSGA2: $(OBJS)/individu.o $(OBJS)/NSGA-2.o $(OBJS)/individuFactory.o $(OBJS)/output.o
	$(CPP) $(DBG) -o $@ $^ $(LDFLAGS) $(CFLAGS) $(OPTIMIZATION)

$(OBJS)/individuFactory.o: $(SRCS)/individuFactory.cpp $(HEADS)/individu.hpp
	$(CPP) $(DBG) -o $@ -c $< $(CFLAGS) $(OPTIMIZATION)

$(OBJS)/individu.o: $(SRCS)/individu.cpp
	$(CPP) $(DBG) -o $@ -c $< $(CFLAGS) $(OPTIMIZATION)

$(OBJS)/output.o: $(SRCS)/output.cpp $(HEADS)/individu.hpp
	$(CPP) $(DBG) -o $@ -c $< $(CFLAGS) $(OPTIMIZATION)

$(OBJS)/NSGA-2.o: $(SRCS)/NSGA-2.cpp $(HEADS)/individu.hpp $(HEADS)/individuFactory.hpp $(HEADS)/output.hpp
	$(CPP) $(DBG) -o $@ -c $< $(CFLAGS) $(OPTIMIZATION)

clean:
	rm -rf $(OBJS)/*.o

mrproper: clean
	rm -rf $(AIMS)

blank:
	tput clear
