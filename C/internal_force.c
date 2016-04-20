#include <math.h>
#include <petsc.h>
#include "global_variables.h"
#include "internal_force.h"
#include "helpers.h"
#include "shape.h"
#include "integration.h"
#include "material.h"


int InternalForce(Vec* dofs, Vec* fglo) {
	PetscErrorCode ierr;
	int i, j, l, ele, row, intpt;
	double kappa = materialprops[1];
	// Initialize
	ierr = VecZeroEntries(*fglo); CHKERRQ(ierr);

	int npt = IntNum(nsd, nen, 0); // Number of integration points
	double xilist[npt][nsd]; // Natural coordinates of the integration points
	IntPoints(nsd, nen, npt, xilist);
	double weights[npt]; // Weights of the integration points
	IntWeights(nsd, nen, npt, weights);

	double elecoord[nen][nsd]; // Coordinates of the nodes in the element
	double eledof[nen][nsd]; // Displacements of the nodes in the element
	double pressure; // Pressure of the element
	double fele[nen*nsd+1]; // Element internal force
	
	// Loop over elements
	for (ele = 0; ele < nel; ++ele) {
		// Extract coords of coords of the nodes and the dofs of the element
		for (i = 0; i < nen; ++i) {
			for (j = 0; j < nsd; ++j) {
				elecoord[i][j] = coords[connect[ele][i]-1][j];
				int id = nsd*(connect[ele][i]-1)+j;
				ierr = VecGetValues(*dofs, 1, &id, &eledof[i][j]); CHKERRQ(ierr);
				//eledof[i][j] = dofs[nsd*(connect[ele][i]-1)+j];
			}
		}
		int id = nsd*nn + ele;
		ierr = VecGetValues(*dofs, 1, &id, &pressure); CHKERRQ(ierr);
		//pressure = dofs[nsd*nn+ele];
		for (i = 0; i < nen*nsd+1; ++i) fele[i] = 0.;
		// Loop over the integration points
		for (intpt = 0; intpt < npt; ++intpt) {
			double N[nen]; // The nodal values of the shape functions
			ShapeFun(nsd, xilist[intpt], nen, N);
			double intcoord[nsd];
			MatMultTranspose_(nen, nsd, elecoord, N, intcoord); // elecoord^T * N
			double dNdxi[nen][nsd];
			ShapeDer(nsd, xilist[intpt], nen, dNdxi);
			double dxdxi[nsd][nsd];
			MatTransposeMatMult_(nsd, nen, nsd, elecoord, dNdxi, dxdxi); // elecoord^T * dNdxi
			double det = Determinant(nsd, dxdxi);
			double dxidx[nsd][nsd];
			InvertDirect(nsd, dxdxi, dxidx);
			double dNdx[nen][nsd];
			MatMatMult_(nen, nsd, nsd, dNdxi, dxidx, dNdx); // dNdxi * dxidx
			double F[nsd][nsd];
			MatTransposeMatMult_(nsd, nen, nsd, eledof, dNdx, F); // I + eledof^T * dNdx
			MatMatAdd_(nsd, nsd, F, i3, F);
			double J = Determinant(nsd, F);
			double stress[nsd][nsd];
			KirchhoffStress(intcoord, F, pressure, stress);
			double Finv[nsd][nsd];
			InvertDirect(nsd, F, Finv);
			double dNdy[nen][nsd];
			MatMatMult_(nen, nsd, nsd, dNdx, Finv, dNdy); // dNdx * Finv
			/*
			printf("intcoord:\n");
			for (i = 0; i < nsd; ++i) printf("%12.4e ", intcoord[i]);
			printf("\n");
			printf("F:\n");
			for (i = 0; i < nsd; ++i) {
				for (j = 0; j < nsd; ++j) {
					printf("%12.4e ", F[i][j]);
				}
				printf("\n");
			}
			printf("pressure: %12.4e\n", pressure);
			*/
			
			// Element internal force
			for (l = 0; l < nen; ++l) {
				for (i = 0; i < nsd; ++i) {
					row = nsd*l + i;
					for (j = 0; j < nsd; ++j) {
						fele[row] += stress[i][j]*dNdy[l][j]*weights[intpt]*det;
					}
				}
			}
			fele[nen*nsd] += (J - 1 - pressure/kappa)*weights[intpt]*det;
		}
		// Scatter
		for (l = 0; l < nen; ++l) {
			for (i = 0; i < nsd; ++i) {
				row = nsd*(connect[ele][l] - 1) + i;
				ierr = VecSetValue(*fglo, row, fele[nsd*l+i], ADD_VALUES); CHKERRQ(ierr);
				//fglo[row] += fele[nsd*l + i];
			}
		}
		row = nn*nsd + ele;
		ierr = VecSetValue(*fglo, row, fele[nsd*nen], ADD_VALUES); CHKERRQ(ierr);
		//fglo[row] += fele[nsd*nen];
	}
	ierr = VecAssemblyBegin(*fglo); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*fglo); CHKERRQ(ierr);
	return 0;
}