#include <math.h>
#include "global_variables.h"
#include "material.h"
#include "helpers.h"

// Computes the Kirchhoff stress at the integration point
void KirchhoffStress(double intcoord[nsd], double F[nsd][nsd], double pressure, double stress[nsd][nsd]) {
	int i, j;
	for (i = 0; i < nsd; ++i) {
		for (j = 0; j < nsd; ++j) {
			stress[i][j] = 0.0;
		}
	}
	// I1 I2 are actually I1bar I2bar respectively
	double I1, I2, J;
	double Bbar[nsd][nsd], BB[nsd][nsd];
	J = Determinant(nsd, F);
	MatMatTransposeMult_(nsd, nsd, nsd, F, F, Bbar);
	MatScale_(pow(J, -2/3.), nsd, nsd, Bbar);
	MatMatMult_(nsd, nsd, nsd, Bbar, Bbar, BB);

	if (nsd == 2) I1 = pow(J, -2/3.);
	else if (nsd == 3) I1 = 0.;
	for (i = 0; i < nsd; ++i) I1 += Bbar[i][i];
	I2 = pow(I1, 2.);
	for (i = 0; i < nsd; ++i) {
		for (j = 0; j < nsd; ++j) {
			I2 -= pow(Bbar[i][j], 2.);
		}
	}
	if (nsd == 2) I2 -= pow(J, -4/3.);
	I2 /= 2.;

	if (materialtype == 1) {
		// Mooney-Rivlin Model
		//double kappa = materialprops[1];
		double mu1 = materialprops[2];
		double mu2 = materialprops[3];
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				stress[i][j] = -(mu1*I1 + 2*mu2*I2)*DELTA(i,j)/3. - mu2*BB[i][j] + (mu1+mu2*I1)*Bbar[i][j] + J*pressure*DELTA(i,j);
			}
		}
	} else if (materialtype == 2) {
		// Yeoh Model
		//double kappa = materialprops[1];
		double c1 = materialprops[2];
		double c2 = materialprops[3];
		double c3 = materialprops[4];
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				stress[i][j] = (2*c1 + 4*c2*(I1 - 3) + 6*c3*(I1-3)*(I1-3))*(Bbar[i][j] - 1/3.*I1*DELTA(i, j)) + J*pressure*DELTA(i, j);
			}
		}
	} else if (materialtype == 3) {
		// Holzapfel-Gasser-Ogden Model
		//double kappa = materialprops[1];
		double mu1 = materialprops[2];
		double k1 = materialprops[3];
		double k2 = materialprops[4];

		// These material constants are hard-coded
		double R = sqrt( pow(intcoord[0], 2) + pow(intcoord[1], 2) );
		double beta;
		if (R <= 0.97e-3) {
			beta = 29*3.14159265/180.0;
			k1 = 2.3632e3;
			k2 = 0.8393;
			mu1 = 3000;
		} else {
			beta = 62*3.14159265/180.0;
			k1 = 0.5620e3;
			k2 = 0.7112;
			mu1 = 300;
		}

		double a0[nsd]; // assign manually
		double g0[nsd]; // assign manually
		if (nsd == 3) {
			a0[0] =  cos(beta)*intcoord[1]/R;
			a0[1] = -cos(beta)*intcoord[0]/R;
			a0[2] =  sin(beta);
			g0[0] =  a0[0];
			g0[1] =  a0[1];
			g0[2] = -a0[2];
		} else if (nsd == 2) {
			a0[0] =  cos(beta)*intcoord[1]/R;
			a0[1] = -cos(beta)*intcoord[0]/R;
			g0[0] =  a0[0];
			g0[1] =  a0[1];
		}

		double a[nsd], g[nsd], temp[nsd];
		double C[nsd][nsd];
		double I4, I6, lambda4, lambda6;
		MatTransposeMatMult_(nsd, nsd, nsd, F, F, C);

		MatMult_(nsd, nsd, C, a0, temp);
		I4 = DotProduct(nsd, a0, temp);
		if (nsd == 2) I4 += sin(beta)*sin(beta);
		lambda4 = sqrt(I4);
		MatMult_(nsd, nsd, F, a0, a);
		VecScale_(1.0/lambda4, nsd, a);
		I4 *= pow(J, -2/3.);

		MatMult_(nsd, nsd, C, g0, temp);
		I6 = DotProduct(nsd, g0, temp);
		if (nsd == 2) I6 += sin(beta)*sin(beta);
		lambda6 = sqrt(I6);
		MatMult_(nsd, nsd, F, g0, g);
		VecScale_(1.0/lambda6, nsd, g);
		I6 *= pow(J, -2/3.);

		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				stress[i][j] = -1/3.*mu1*I1*DELTA(i,j) + mu1*Bbar[i][j] + J*pressure*DELTA(i,j);
				// Anisotropic Part
				if (lambda4 > 1.0) {
					stress[i][j] += 2*k1*( (I4-1)*exp(k2*(I4-1)*(I4-1))*(a[i]*a[j] - 1/3.*DELTA(i,j))*I4 );
				}
				if (lambda6 > 1.0) {
					stress[i][j] += 2*k1*( (I6-1)*exp(k2*(I6-1)*(I6-1))*(g[i]*g[j] - 1/3.*DELTA(i,j))*I6 );
				}
			}
		}
	}
}

// Computes the 4th-order spatial tensor of elasiticity at the integration point
void MaterialStiffness(double intcoord[nsd], double F[nsd][nsd], double pressure, double mstiff[nsd][nsd][nsd][nsd]) {
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

	// I1 I2 are actually I1bar I2bar respectively
	double I1, I2, J;
	double Bbar[nsd][nsd], BB[nsd][nsd];
	J = Determinant(nsd, F);
	MatMatTransposeMult_(nsd, nsd, nsd, F, F, Bbar);
	MatScale_(pow(J, -2/3.), nsd, nsd, Bbar);
	MatMatMult_(nsd, nsd, nsd, Bbar, Bbar, BB);

	if (nsd == 2) I1 = pow(J, -2/3.);
	else if (nsd == 3) I1 = 0.;
	for (i = 0; i < nsd; ++i) I1 += Bbar[i][i];
	I2 = pow(I1, 2.);
	for (i = 0; i < nsd; ++i) {
		for (j = 0; j < nsd; ++j) {
			I2 -= pow(Bbar[i][j], 2);
		}
	}
	if (nsd == 2) I2 -= pow(J, -4/3.);
	I2 /= 2.;

	if (materialtype == 1) {
		// Mooney-Rivlin Model
		//double kappa = materialprops[1];
		double mu1 = materialprops[2];
		double mu2 = materialprops[3];
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				for (k = 0; k < nsd; ++k) {
					for (l = 0; l < nsd; ++l) {
						mstiff[i][j][k][l] = mu2*(2*Bbar[i][j]*Bbar[k][l] - (Bbar[i][k]*Bbar[j][l] + Bbar[i][l]*Bbar[j][k])) \
						- 2/3.*(mu1 + 2*mu2*I1)*(Bbar[k][l]*DELTA(i, j) + Bbar[i][j]*DELTA(k, l)) \
						+ 4/3.*mu2*(BB[i][j]*DELTA(k, l) + BB[k][l]*DELTA(i, j)) \
						+ 2/9.*(mu1*I1 + 4*mu2*I2)*DELTA(i, j)*DELTA(k, l) \
						+ 1/3.*(mu1*I1 + 2*mu2*I2)*(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)) \
						+ (DELTA(i, j)*DELTA(k, l) - (DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)))*J*pressure;
					}
				}
			}
		}
	} else if (materialtype == 2) {
		// Yeoh Model
		//double kappa = materialprops[1];
		double c1 = materialprops[2];
		double c2 = materialprops[3];
		double c3 = materialprops[4];
		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				for (k = 0; k < nsd; ++k) {
					for (l = 0; l < nsd; ++l) {
						mstiff[i][j][k][l] = (8*c2 + 24*c3*(I1 - 3))*Bbar[i][j]*Bbar[k][l] \
						- (4/3.*c1 + (16/3.*I1 - 8)*c2 + 12*(I1 - 1)*(I1 - 3)*c3)*(DELTA(i, j)*Bbar[k][l] + DELTA(k, l)*Bbar[i][j]) \
						+ (4/9.*I1*c1 + (16/9.*I1*I1 - 8/3.*I1)*c2 + 4*I1*(I1 - 1)*(I1 - 3)*c3)*DELTA(i, j)*DELTA(k, l) \
						+ (2/3.*I1*c1 + 4/3.*I1*(I1 - 3)*c2 + 2*I1*(I1 - 3)*(I1 - 3)*c3)* \
							(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)) \
						+ (DELTA(i, j)*DELTA(k, l) - (DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)))*J*pressure;
					}
				}
			}
		}
	} else if (materialtype == 3) {
		// Holzapfel-Gasser-Ogden Model
		//double kappa = materialprops[1];
		double mu1 = materialprops[2];
		double k1 = materialprops[3];
		double k2 = materialprops[4];

		// These material constants are hard-coded
		double R = sqrt( pow(intcoord[0], 2) + pow(intcoord[1], 2) );
		double beta;
		if (R <= 0.97e-3) {
			beta = 29*3.14159265/180.0;
			k1 = 2.3632e3;
			k2 = 0.8393;
			mu1 = 3000;
		} else {
			beta = 62*3.14159265/180.0;
			k1 = 0.5620e3;
			k2 = 0.7112;
			mu1 = 300;
		}

		double a0[nsd]; // assign manually
		double g0[nsd]; // assign manually
		if (nsd == 3) {
			a0[0] =  cos(beta)*intcoord[1]/R;
			a0[1] = -cos(beta)*intcoord[0]/R;
			a0[2] =  sin(beta);
			g0[0] =  a0[0];
			g0[1] =  a0[1];
			g0[2] = -a0[2];
		} else if (nsd == 2) {
			a0[0] =  cos(beta)*intcoord[1]/R;
			a0[1] = -cos(beta)*intcoord[0]/R;
			g0[0] =  a0[0];
			g0[1] =  a0[1];
		}

		double a[nsd], g[nsd], temp[nsd];
		double C[nsd][nsd];
		double I4, I6, lambda4, lambda6;

		MatTransposeMatMult_(nsd, nsd, nsd, F, F, C);
		MatMult_(nsd, nsd, C, a0, temp);
		I4 = DotProduct(nsd, a0, temp);
		if (nsd == 2) I4 += sin(beta)*sin(beta);
		lambda4 = sqrt(I4);
		MatMult_(nsd, nsd, F, a0, a);
		VecScale_(1.0/lambda4, nsd, a);
		I4 *= pow(J, -2/3.);

		MatMult_(nsd, nsd, C, g0, temp);
		I6 = DotProduct(nsd, g0, temp);
		if (nsd == 2) I6 += sin(beta)*sin(beta);
		lambda6 = sqrt(I6);
		MatMult_(nsd, nsd, F, g0, g);
		VecScale_(1.0/lambda6, nsd, g);
		I6 *= pow(J, -2/3.);

		double der[2][2]; 
		// 1st row: 1st order derivative; 2nd row: 2nd order derivative; 1st col: w.r.t. I4; 2nd col: w.r.t. I6
		der[0][0] = k1 * (I4 - 1.) * exp(k2*(I4 - 1.)*(I4 - 1.));
		der[0][1] = k1 * (I6 - 1.) * exp(k2*(I6 - 1.)*(I6 - 1.));
		der[1][0] = k1 * (1. + 2*k2*(I4 - 1.)*(I4 - 1.))*exp(k2*(I4-1.)*(I4-1.));
		der[1][1] = k1 * (1. + 2*k2*(I6 - 1.)*(I6 - 1.))*exp(k2*(I6-1.)*(I6-1.));

		for (i = 0; i < nsd; ++i) {
			for (j = 0; j < nsd; ++j) {
				for (k = 0; k < nsd; ++k) {
					for (l = 0; l < nsd; ++l) {
						mstiff[i][j][k][l] = - 2/3.*mu1*(Bbar[k][l]*DELTA(i, j) + Bbar[i][j]*DELTA(k, l)) \
						+ 2/9.*mu1*I1*DELTA(i, j)*DELTA(k, l) \
						+ 1/3.*mu1*I1*(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)) \
						+ (DELTA(i, j)*DELTA(k, l) - (DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k)))*J*pressure;
						// Anisotropic part 1
						if (lambda4 > 1.0) {
							mstiff[i][j][k][l] += 4*I4*I4*der[1][0]*a[i]*a[j]*a[k]*a[l] \
							- 4/3.*(I4*der[1][0] + der[0][0])*I4*(a[i]*a[j]*DELTA(k, l) + a[k]*a[l]*DELTA(i, j)) \
							+ (4/9.*I4*I4*der[1][0] + 4/9.*I4*der[0][0])*DELTA(i, j)*DELTA(k, l) \
							+ 2/3.*I4*der[0][0]*(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k));
						}
						// Anisotropic part 2
						if (lambda6 > 1.0) {
							mstiff[i][j][k][l] += 4*I6*I6*der[1][1]*g[i]*g[j]*g[k]*g[l] \
							- 4/3.*(I6*der[1][1] + der[0][1])*I6*(g[i]*g[j]*DELTA(k, l) + g[k]*g[l]*DELTA(i, j)) \
							+ (4/9.*I6*I6*der[1][1] + 4/9.*I6*der[0][1])*DELTA(i, j)*DELTA(k, l) \
							+ 2/3.*I6*der[0][1]*(DELTA(i, k)*DELTA(j, l) + DELTA(i, l)*DELTA(j, k));
						}
					}
				}
			}
		}
	}
}

