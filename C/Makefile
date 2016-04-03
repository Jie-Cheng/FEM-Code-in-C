LIBS = /usr/local/lib/liblapack.a  /usr/local/lib/libblas.a /usr/local/lib/libmetis.a /usr/local/lib/libma57.a
OBJECTS = read.o helpers.o face.o integration.o shape.o material.o internal_force.o \
		  external_force.o tangent_stiffness.o symmetric_solver.o main.o
FC = gfortran
CC = gcc
CPP = g++

main: $(OBJECTS)
	$(FC) -o main $(OBJECTS) $(LIBS)

read.o: read.cpp
	$(CPP) -c read.cpp

symmetric_solver.o: symmetric_solver.f90
	$(FC) -c symmetric_solver.f90

%.o: %.c
	$(CC) -c $<

clean:
	rm -f *.o main
