INCLUDE = -I/usr/include/ 
#-I/usr/include/python2.7/
LIBDIR  = -L/usr/lib

FLAGS = -Wall
CC = g++                                  # change to gcc if using C
CFLAGS = $(FLAGS) $(INCLUDE)
# LIBS = -lpython2.7


#All: clean eigenSolver.o visual.o EigenProj                             # change your_app.

Solver: main.o
		$(CC) $(CFLAGS) -o $@ $(LIBDIR) $^ $(LIBS) # The initial white space is a tab

# main.o: main.cpp Solver.h Solver1d.h CallPy.h spectrum.h function.h
main.o: main.cpp Solver.h Solver2d.h spectrum.h function.h
		$(CC) -c main.cpp $(LIBS) $(INCLUDE)

clean:  
	rm -rf *.o  
	rm -rf Solver

