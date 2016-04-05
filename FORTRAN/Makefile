LIBS = /usr/local/lib/liblapack.a  /usr/local/lib/libblas.a /usr/local/lib/libmetis.a /usr/local/lib/libma57.a
OBJECTS = read_file.o integration.o face.o shapefunction.o material.o externalforce.o output.o \
	      symmetric_solver.o mass.o \
		  tangentstiffness.o internalforce.o solidsolver.o
MODULES = read_file.mod integration.mod face.mod shapefunction.mod material.mod externalforce.mod \
          symmetric_solver.mod mass.mod \
		  tangentstiffness.mod internalforce.mod output.mod
FC = gfortran

main: $(OBJECTS)
	$(FC) -o main $(OBJECTS) $(LIBS)
			
%.mod: %.f90
	$(FC) -c $<
	
%.o: %.f90
	$(FC) -c $<
	
clean:
	rm -f *.o *.mod main
