#include <stdlib.h>
#include <math.h>
#include "material.h"
#include "helpers.h"

void MaterialStiffness(int nsd, double intcoord[nsd], double F[nsd][nsd], \
	double pressure, int materialtype, double* materialprops, double mstiff[nsd][nsd][nsd][nsd]) {
	int i, j, k, l;

	for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				for (k = 0; k < nsd; ++k) {
					for (l = 0; l < nsd; ++l) {
						mstiff[i][j][k][l] = 0.0;
					}
				}
			}
		}
	}

	// I1 I2 are actually I1bar I2bar respectively
	double I1, I2, J;
	double B[nsd][nsd], Bbar[nsd][nsd], BB[nsd][nsd];

	J = Determinant(F, nsd);
	MatMul('N','T', nsd, nsd, nsd, 1.0, 0.0, F[0], F[0], B[0]);
	MatMul('N','T', nsd, nsd, nsd, pow(J, -2/3.), 0.0, F[0], F[0], B[0]);
	MatMul('N','N', nsd, nsd, nsd, 1.0, 0.0, Bbar[0], Bbar[0], BB[0]);

	if (nsd == 2) I1 = 1.;
	else if (nsd == 3) I1 = 0.;
	for (i = 0; i < nsd; ++i) I1 += Bbar[i][i];
		I2 = pow(I1, 2);
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				I2 -= pow(Bbar[i][j], 2)
			}
		}
	}
	if (nsd == 2) I2 -= pow(J, -4/3.);
	I2 /= 2.;

	if (materialtype == 1) {
		// Mooney-Rivlin Model
		double mu1 = materialprops[2];
		double mu2 = materialprops[3];
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				for (k = 0; k < nsd; ++k) {
					for (l = 0; l < nsd; ++l) {
						mstiff[i][j][k][l] += mu2*(2*Bbar[i][j]*Bbar[k][l] - (Bbar[i][k]*Bbar[j][l] + Bbar[i][l]*Bbar[j][k])) \
						- 2/3.*(mu1 + 2*mu2*I1)*(Bbar[k][l]*DELTA(i, j) + Bbar[i][j]*DELTA(k, l)) \
						+ 4/3.*mu2*(BB[i][j]*DELTA(k, l) + BB[k][l]*DELTA(i, j)) \
						+ 2/9.*(mu1*I1 + 4*mu2*I2)*(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)) \
						+ (DELTA(i, j)*DELTA(k, l) - (DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)))*J*pressure;
					}
				}
			}
		}
	} else if (materialtype == 2) {
		// Yeoh Model
		double c1 = materialprops[2];
		double c2 = materialprops[3];
		double c3 = materialprops[4];
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				for (k = 0; k < nsd; ++k) {
					for (l = 0; l < nsd; ++l) {
						mstiff[i][j][k][l] += (8*c2 + 24*c3*(I1 - 3))*Bbar[i][j]*Bbar[k][l] \
						- (4/3.*c1 + (16/3.*I1 - 8)*c2 + 12*(I1 - 1)*(I1 - 3)*c3)*(DELTA(i, j)*Bbar[k][l] + DELTA(k, l)*Bbar[i][j]) \
						+ (4/9.*I1*c1 + (16/9.*I1*I1 - 8/3.*I1)*c2 + 4*I1*(I1 - 3)*(I1 - 3)*c3)*DELTA(i, j)*DELTA(k, l) \
						+ (2/3.*I1*c1 + 4/3.*I1*(I1 - 3)*c2 + (2*I1*(I1 - 3)*(I1 - 3))*c3)* \
							(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k))
						+ (DELTA(i, j)*DELTA(k, l) - (DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)))*J*pressure;
					}
				}
			}
		}
	} else if (materialtype == 3) {
		// Holzapfel-Gasser-Ogden Model
		// These material constants are hard-coded
		double beta = 50*M_PI/180.0;
		double R = sqrt( (1 + pow(tan(beta), 2)) * (pow(intcoord[0], 2) + pow(intcoord[1], 2)) );
		double a0[nsd]; // assign manually
		double g0[nsd]; // assign manually
		double k1 = 0.5158e6;
		double k2 = 0.9;

		double a[nsd], g[nsd];
		double C[nsd][nsd];
		double I4, I6, lambda;

		MatMul('T','N', nsd, nsd, nsd, 1.0, 0.0, F[0], F[0], B[0]);
		


	}

}