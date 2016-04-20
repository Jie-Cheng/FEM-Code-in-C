#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <petscksp.h>
#include "global_variables.h"
#include "read.h"
#include "internal_force.h"
#include "external_force.h"
#include "tangent_stiffness.h"

int Statics();
int Debug();

int main(int argc,char* argv[]) {
	// Initialization
	ReadInput();
	ReadMesh();
	Echo();
	CountNonzerosInSerial();

	PetscErrorCode ierr;
	PetscMPIInt mycomm;

	PetscInitialize(&argc, &argv, (char*)0, (char*)0);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &mycomm); CHKERRQ(ierr);
	if (mycomm != 1) 
		SETERRQ(PETSC_COMM_WORLD, 1, "For now this program can only run in serial!");
	Statics();
	//Debug();
	

	ierr = PetscFinalize();

	int i;
	for (i = 0; i < nn; ++i) free(coords[i]);
	free(coords);
	for (i = 0; i < nel; ++i) free(connect[i]);
	free(connect);
	for (i = 0; i < bc_size; ++i) free(bc_num[i]);
	free(bc_num);
	free(bc_val);
	for (i = 0; i < load_size; ++i) {
		free(load_num[i]);
		free(load_val[i]);
	}
	free(load_num);
	free(load_val);
	free(share);
	
    return 0;
}

int Statics() {
	int ndofs = nn*nsd + nel;
	
	PetscErrorCode ierr;
	Vec Fint, Fext, R, w, w1, dw;
	Mat A;
	ierr = VecCreate(PETSC_COMM_WORLD, &Fint); CHKERRQ(ierr);
	ierr = VecSetSizes(Fint, PETSC_DECIDE, ndofs); CHKERRQ(ierr);
	ierr = VecSetFromOptions(Fint); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &Fext); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &R); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &w); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &w1); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &dw); CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, ndofs, ndofs, PETSC_DEFAULT, nnz, &A); CHKERRQ(ierr);
	//ierr = MatSetType(A, MATSEQSBAIJ); CHKERRQ(ierr);
	//ierr = MatSetUp(A);CHKERRQ(ierr);
	KSP ksp; // Linear solver context
	PC pc; // Preconditioner context
	//PetscReal tolorance = 1.0e-14;
	double constraint = 0.0;
	int i, row, gmres_nit;

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = VecSet(w, 0.0); CHKERRQ(ierr);
	ierr = VecSet(w1, 0.0); CHKERRQ(ierr);
	ierr = VecSet(dw, 0.0); CHKERRQ(ierr);
	ierr = VecSet(Fint, 0.0); CHKERRQ(ierr);
	ierr = VecSet(Fext, 0.0); CHKERRQ(ierr);
	ierr = VecSet(R, 0.0); CHKERRQ(ierr);

	nprint = 1;
	nsteps = 1;
	dt = 1.;
	step = 0;
	double load_factor = 0.0;
	double increment = inistep;
	
	// output

	if (load_type != 1) ExternalTraction(&Fext); // Traction

	while (load_factor < 1.0) {
		step++;
		if (load_factor + increment > 1.0) increment = 1.0 - load_factor;
		load_factor += increment;
		ierr = VecCopy(w1, w); CHKERRQ(ierr);
		double err1 = 1.0;
		double err2 = 1.0;
		int nit = 0;
		printf("Step = %4d, LoadFactor = %12.4e\n", step, load_factor);
		while ((err1 > tol || err2 > tol) && nit < maxit) {
			nit++;
			InternalForce(&w, &Fint);
			if (load_type == 1) ExternalPressure(&w, &Fext); // pressure load
			InternalTangent(&w, &A);
			//ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
			// R =  Fint - load_factor*Fext
			ierr = VecCopy(Fext, R); CHKERRQ(ierr);
			ierr = VecAYPX(R, -load_factor, Fint); CHKERRQ(ierr);
			for (i = 0; i < bc_size; ++i) {
				row = nsd*(bc_num[i][0] - 1) + bc_num[i][1] - 1;
				ierr = VecGetValues(w, 1, &row, &constraint); CHKERRQ(ierr);
				constraint -= bc_val[i];
				ierr = MatSetValue(A, row, row, penalty, ADD_VALUES); CHKERRQ(ierr);
				ierr = VecSetValue(R, row, penalty*constraint, ADD_VALUES); CHKERRQ(ierr);
			}
			ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
			ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
			ierr = VecAssemblyBegin(R); CHKERRQ(ierr);
			ierr = VecAssemblyEnd(R); CHKERRQ(ierr);
			ierr = VecScale(R, -1.0); CHKERRQ(ierr);

			ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
			ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
			ierr = PCSetType(pc, PCKSP); CHKERRQ(ierr);
			ierr = KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
			ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
			printf("Begin to solve\n");
			ierr = KSPSolve(ksp, R, dw); CHKERRQ(ierr);
			printf("End of solve\n");
			KSPConvergedReason reason;
			ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
			printf("REASON: %d\n", reason);

			ierr = KSPGetIterationNumber(ksp, &gmres_nit); CHKERRQ(ierr);
			ierr = VecAYPX(w, 1.0, dw); CHKERRQ(ierr);

			double temp;
			//ierr = VecView(w, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
			ierr = VecDot(dw, dw, &err1); CHKERRQ(ierr);
			ierr = VecDot(w, w, &temp); CHKERRQ(ierr);
			err1 = sqrt(err1/temp);
			ierr = VecDot(R, R, &err2); CHKERRQ(ierr);
			err2 = sqrt(err2)/ndofs;
			ierr = VecZeroEntries(dw); CHKERRQ(ierr);
			
			printf("Iteration = %4d     GMRES_Iteration = %4d     Err1 = %12.4e,     Err2 = %12.4e\n",\
				nit, gmres_nit, err1, err2);

		}
		if (nit == maxit) {
			ierr = VecZeroEntries(w); CHKERRQ(ierr);
			ierr = VecCopy(w, w1); CHKERRQ(ierr);
			load_factor = load_factor - increment;
			increment /= 2.;
		} else if (nit > 6) {
			increment *= (1.0 - adjust);
		} else if (nit < 6) {
			increment *= (1.0 + adjust);
		}
	}
	step = nprint;

	ierr = VecView(w, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

	ierr = VecDestroy(&Fint); CHKERRQ(ierr);
	ierr = VecDestroy(&Fext); CHKERRQ(ierr);
	ierr = VecDestroy(&R); CHKERRQ(ierr);
	ierr = VecDestroy(&w); CHKERRQ(ierr);
	ierr = VecDestroy(&w1); CHKERRQ(ierr);
	ierr = VecDestroy(&dw); CHKERRQ(ierr);
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

	// output
	return 0;
}


int Debug() {
	int ndofs = nn*nsd + nel;
	
	PetscErrorCode ierr;
	Vec Fint, Fext, R, w, w1, dw;
	Mat A;
	ierr = VecCreate(PETSC_COMM_WORLD, &Fint); CHKERRQ(ierr);
	ierr = VecSetSizes(Fint, PETSC_DECIDE, ndofs); CHKERRQ(ierr);
	ierr = VecSetFromOptions(Fint); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &Fext); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &R); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &w); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &w1); CHKERRQ(ierr);
	ierr = VecDuplicate(Fint, &dw); CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, ndofs, ndofs, PETSC_DEFAULT, nnz, &A); CHKERRQ(ierr);

	int i;
	int ix[ndofs];
	double val[ndofs];
	for (i = 0; i < ndofs; ++i) {
		ix[i] = i;
		val[i] = 0.01*i;
	}

	ierr = VecSetValues(w, ndofs, ix, val, INSERT_VALUES); CHKERRQ(ierr);
	InternalForce(&w, &Fint);
	InternalTangent(&w, &A);
	ExternalPressure(&w, &Fext);

	ierr = VecView(Fint, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	ierr = VecView(Fext, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


	ierr = VecDestroy(&Fint); CHKERRQ(ierr);
	ierr = VecDestroy(&Fext); CHKERRQ(ierr);
	ierr = VecDestroy(&R); CHKERRQ(ierr);
	ierr = VecDestroy(&w); CHKERRQ(ierr);
	ierr = VecDestroy(&w1); CHKERRQ(ierr);
	ierr = VecDestroy(&dw); CHKERRQ(ierr);
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	return 0;
}
