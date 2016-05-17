OBJECTS = global_variables.o read.o shape.o face.o integration.o helpers.o \
		  external_force.o internal_force.o tangent_stiffness.o material.o \
		  analyze_nonzero_pattern.o main.o
CC = gcc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: $(OBJECTS)
	-${CLINKER} -o main $^ ${PETSC_KSP_LIB}

analyze_nonzero_pattern.o: analyze_nonzero_pattern.cpp
	g++ -c -g $<

%.o: %.c
	$(PETSC_COMPILE) -c -g $<

# $< = first input
# $^ = inputs
# $@ = outputs