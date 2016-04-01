#include <math.h>
#include "internal_force.h"
#include "helpers.h"
#include "shape.h"
#include "integration.h"
#include "material.h"

void InternalForce(int nsd, int nn, int nel, int nen, double dofs[nn*nsd+nel], \
	double coords[nn][nsd], int connect[nel][nen], int materialtype, \
	double materialprops[5], double fglo[nn*nsd+nel]) {

	int i, j, k, l, ele, row;
	double kappa = materialprops[1];

	for (i = 0; i < nn*nsd+nel; ++i) fglo[i] = 0.;

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
				eledof[i][j] = dofs[nsd*(connect[ele][i]-1)+j];
			}
		}
		pressure = dofs[nsd*nn+ele];
		for (i = 0; i < nen*nsd+1; ++i) fele[i] = 0.;
		// Loop over the integration points
		int intpt;
		for (intpt = 0; intpt < npt; ++intpt) {
			double N[nen]; // The nodal values of the shape functions
			ShapeFun(nen, N, nsd, xilist[intpt]);
			double intcoord[nsd];
			MatMul('T', 'N', nsd, 1, nen, 1.0, 0.0, elecoord[0], N, intcoord); // MatMulVec(elecoord^T, N)
			double dNdxi[nen][nsd];
			ShapeDer(nsd, xilist[intpt], nen, dNdxi);
			double dxdxi[nsd][nsd];
			MatMul('T', 'N', nsd, nsd, nen, 1.0, 0.0, elecoord[0], dNdxi[0], dxdxi[0]); // MatMulMat(elecoord^T, dNdxi)
			double det = Determinant(nsd, dxdxi);
			double dxidx[nsd][nsd];
			for (i = 0; i < nsd; ++i) {
				for (j = 0; j < nsd; ++j) {
					dxidx[i][j] = dxdxi[i][j];
				}
			}
			Inverse(nsd, dxidx[0]);
			double dNdx[nen][nsd];
			MatMul('N', 'N', nen, nsd, nsd, 1.0, 0.0, dNdxi[0], dxidx[0], dNdx[0]); // MatMulMat(dNdxi, dxidx)
			double F[nsd][nsd];
			MatMul('T', 'N', nsd, nsd, nen, 1.0, 0.0, eledof[0], dNdx[0], F[0]); // eye + MatMul(eledof^T, dNdx)
			for (i = 0; i < nsd; ++i) {
				for (j = 0; j < nsd; ++j) {
					F[i][j] += 1.0;
				}
			}
			double J = Determinant(nsd, F);
			double stress[nsd][nsd];
			KirchhoffStress(nsd, intcoord, F, pressure, materialtype, materialprops, stress);
			double Finv[nsd][nsd];
			for (i = 0; i < nsd; ++i) {
				for (j = 0; j < nsd; ++j) {
					Finv[i][j] = F[i][j];
				}
			}
			Inverse(nsd, Finv[0]);
			double dNdy[nen][nsd]; 
			MatMul('N', 'N', nen, nsd, nsd, 1.0, 0.0, dNdx[0], Finv[0], dNdy[0]); // MatMul(dNdx, Finv)
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
				fglo[row] += fele[nsd*l + i];
			}
		}
		row = nn*nsd + ele;
		fglo[row] += fele[nsd*nen];	
	}	
}