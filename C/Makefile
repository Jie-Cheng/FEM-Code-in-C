OBJECTS = global_variables.o read.o shape.o face.o integration.o helpers.o \
		  external_force.o internal_force.o tangent_stiffness.o material.o main.o
CC = gcc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: $(OBJECTS)
	-${CLINKER} -o main $^ ${PETSC_KSP_LIB}

%.o: %.c
	$(PETSC_COMPILE) -c -g $<

myclean:
	rm -f *.o main

# $< = first input
# $^ = inputs
# $@ = outputs