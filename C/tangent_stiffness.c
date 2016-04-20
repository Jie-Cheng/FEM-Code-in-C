#include <math.h>
#include <petscksp.h>
#include "global_variables.h"
#include "tangent_stiffness.h"
#include "helpers.h"
#include "shape.h"
#include "integration.h"
#include "material.h"

int InternalTangent(Vec* dofs, Mat* kglo) {
	PetscErrorCode ierr;

	int i, j, k, l, a, b, ele, row, col, intpt;
	double kappa = materialprops[1];

	ierr = MatZeroEntries(*kglo); CHKERRQ(ierr);

	int npt = IntNum(nsd, nen, 0); // Number of integration points
	double xilist[npt][nsd]; // Natural coordinates of the integration points
	IntPoints(nsd, nen, npt, xilist);
	double weights[npt]; // Weights of the integration points
	IntWeights(nsd, nen, npt, weights);

	double elecoord[nen][nsd]; // Coordinates of the nodes in the element
	double eledof[nen][nsd]; // Displacements of the nodes in the element
	double pressure; // Pressure of the element
	double kele[nen*nsd+1][nen*nsd+1]; // Element internal tangent stiffness

	for (ele = 0; ele < nel; ++ele) {
		for (i = 0; i < nen; ++i) {
			for (j = 0; j < nsd; ++j) {
				elecoord[i][j] = coords[connect[ele][i]-1][j];
				int id = nsd*(connect[ele][i]-1)+j;
				ierr = VecGetValues(*dofs, 1, &id, &eledof[i][j]); CHKERRQ(ierr);
			}
		}
		int id = nsd*nn + ele;
		ierr = VecGetValues(*dofs, 1, &id, &pressure); CHKERRQ(ierr);
		for (i = 0; i < nen*nsd+1; ++i) {
			for (j = 0; j < nen*nsd+1; ++j) {
				kele[i][j] = 0.0;
			}
		}
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
			double Finv[nsd][nsd];
			InvertDirect(nsd, F, Finv);
			double dNdy[nen][nsd];
			MatMatMult_(nen, nsd, nsd, dNdx, Finv, dNdy); // dNdx * Finv
			double J = Determinant(nsd, F);
			double stress[nsd][nsd];
			KirchhoffStress(intcoord, F, pressure, stress);
			double C[nsd][nsd][nsd][nsd];
			MaterialStiffness(intcoord, F, pressure, C);
			for (a = 0; a < nen; ++a) {
				for (i = 0; i < nsd; ++i) {
					row = a*nsd + i;
					for (b = 0; b < nen; ++b) {
						for (k = 0; k < nsd; ++k) {
							col = b*nsd + k;
							for (j = 0; j < nsd; ++j) {
								for (l = 0; l < nsd; ++l) {
									kele[row][col] += (DELTA(i, k)*stress[j][l] + \
										C[i][j][k][l])*dNdy[a][j]*dNdy[b][l]*weights[intpt]*det;
								}
							}
						}
					}
					col = nen*nsd;
					kele[row][col] += dNdy[a][i]*J*weights[intpt]*det;
					kele[col][row] = kele[row][col];
				}
			}
			kele[nen*nsd][nen*nsd] -= 1/kappa*weights[intpt]*det;
		}
		/*
		for (i = 0; i < 25; ++i) {
			for (j = 0; j < 25; ++j) {
				printf("%lf ", kele[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		*/
		// Scatter
		for (a = 0; a < nen; ++a) {
			for (i = 0; i < nsd; ++i) {
				row = nsd*(connect[ele][a]-1) + i;
				for (b = 0; b < nen; ++b) {
					for (k = 0; k < nsd; ++k) {
						col = nsd*(connect[ele][b]-1) + k;
						ierr = MatSetValue(*kglo, row, col, kele[nsd*a+i][nsd*b+k], ADD_VALUES); CHKERRQ(ierr);
						//kglo[row][col] += kele[nsd*a+i][nsd*b+k];
					}
				}
				col = nsd*nn + ele;
				ierr = MatSetValue(*kglo, row, col, kele[nsd*a+i][nsd*nen], ADD_VALUES); CHKERRQ(ierr);
				ierr = MatSetValue(*kglo, col, row, kele[nsd*a+i][nsd*nen], ADD_VALUES); CHKERRQ(ierr);
				//kglo[row][col] += kele[nsd*a+i][nsd*nen];
				//kglo[col][row]  = kglo[row][col];
			}
		}
		ierr = MatSetValue(*kglo, nsd*nn+ele, nsd*nn+ele, kele[nsd*nen][nsd*nen], ADD_VALUES); CHKERRQ(ierr);
		//kglo[nsd*nn+ele][nsd*nn+ele] += kele[nsd*nen][nsd*nen];
	}
	ierr = MatAssemblyBegin(*kglo, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*kglo, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	return 0;
} 