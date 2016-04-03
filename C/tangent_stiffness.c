#include <math.h>
#include "tangent_stiffness.h"
#include "helpers.h"
#include "shape.h"
#include "integration.h"
#include "material.h"

void InternalTangent(int nsd, int nn, int nel, int nen, double dofs[nn*nsd+nel], \
	double coords[nn][nsd], int connect[nel][nen], int materialtype, \
	double materialprops[5], double kglo[nn*nsd+nel][nn*nsd+nel]) {

	int i, j, k, l, a, b, ele, row, col;
	double kappa = materialprops[1];

	for (i = 0; i < nn*nsd+nel; ++i)
		for (j = 0; j < nn*nsd+nel; ++j)
			kglo[i][j] = 0.0;

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
				eledof[i][j] = dofs[nsd*(connect[ele][i]-1)+j];
			}
		}
		pressure = dofs[nsd*nn+ele];
		for (i = 0; i < nen*nsd+1; ++i) {
			for (j = 0; j < nen*nsd+1; ++j) {
				kele[i][j] = 0.0;
			}
		}
		int intpt;
		for (intpt = 0; intpt < npt; ++intpt) {
			double N[nen]; // The nodal values of the shape functions
			ShapeFun(nsd, xilist[intpt], nen, N);
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
			double Finv[nsd][nsd];
			for (i = 0; i < nsd; ++i) {
				for (j = 0; j < nsd; ++j) {
					Finv[i][j] = F[i][j];
				}
			}
			Inverse(nsd, Finv[0]);
			double dNdy[nen][nsd]; 
			MatMul('N', 'N', nen, nsd, nsd, 1.0, 0.0, dNdx[0], Finv[0], dNdy[0]); // MatMul(dNdx, Finv)
			double J = Determinant(nsd, F);
			double stress[nsd][nsd];
			KirchhoffStress(nsd, intcoord, F, pressure, materialtype, materialprops, stress);
			double C[nsd][nsd][nsd][nsd];
			MaterialStiffness(nsd, intcoord, F, pressure, materialtype, materialprops, C);
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
		// Scatter
		for (a = 0; a < nen; ++a) {
			for (i = 0; i < nsd; ++i) {
				row = nsd*(connect[ele][a]-1) + i;
				for (b = 0; b < nen; ++b) {
					for (k = 0; k < nsd; ++k) {
						col = nsd*(connect[ele][b]-1) + k;
						kglo[row][col] += kele[nsd*a+i][nsd*b+k];
					}
				}
				col = nsd*nn + ele;
				kglo[row][col] += kele[nsd*a+i][nsd*nen];
				kglo[col][row]  = kglo[row][col];
			}
		}
		kglo[nsd*nn+ele][nsd*nn+ele] += kele[nsd*nen][nsd*nen];
	}
} 