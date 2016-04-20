#include "helpers.h"

double i3[3][3] = {
	{1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0}
};

double i2[2][2] = {
	{1.0, 0.0},
	{0.0, 1.0}
};

// Find the cross product of two vectors.
void Cross3d(double a[3], double b[3], double c[3]) {
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

// Determinant of 2X2 or 3X3 matrix
double Determinant(const int n, double a[n][n]) {
	double det = 0.0;
	if (n == 2) {
		det = a[0][0]*a[1][1] - a[0][1]*a[1][0];
	} else if (n == 3) {
		det = a[0][0]*a[1][1]*a[2][2] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] \
		+  a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0];
	}
	return det;
}

// Invert a 2d or 3d square matrix using direct method
void InvertDirect(const int n, double A[n][n], double B[n][n]) {
	double det = Determinant(n, A);
	if (n == 2) {
		B[0][0] =  A[1][1]/det;
		B[0][1] = -A[0][1]/det;
		B[1][0] = -A[1][0]/det;
		B[1][1] =  A[0][0]/det;
	} else if (n == 3) {
		B[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1])/det;
		B[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2])/det;
		B[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])/det;
		B[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2])/det;
		B[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0])/det;
		B[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])/det;
		B[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0])/det;
		B[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1])/det;
		B[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0])/det;
	}
}

// Dot product
double DotProduct(const int n, double a[n], double b[n]) {
	int i;
	double result = 0.0;
	for (i = 0; i < n; ++i) {
		result += a[i]*b[i];
	}
	return result;
}

// v = alpha * v
void VecScale_(const double alpha, const int n, double v[n]) {
	int i;
	for (i = 0; i < n; ++i) {
		v[i] = alpha*v[i];
	}
}
// B = alpha * A
void MatScale_(const double alpha, const int m, const int n, double A[m][n]) {
	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			A[i][j] *= alpha;
		}
	}
}

// b = A * x
void MatMult_(const int m, const int n, double A[m][n], double x[n], double b[m]) {
	int i, j;
	for (i = 0; i < m; ++i) {
		b[i] = 0.0;
		for (j = 0; j < n; ++j) {
			b[i] += A[i][j]*x[j];
		}
	}
}

// b = A^T * x
void MatMultTranspose_(const int m, const int n, double A[m][n], double x[m], double b[n]) {
	int i, j;
	for (i = 0; i < n; ++i) {
		b[i] = 0.0;
		for (j = 0; j < m; ++j) {
			b[i] += A[j][i]*x[j];
		}
	}
}

// C = A * B
void MatMatMult_(const int m, const int k, const int n, \
	double A[m][k], double B[k][n], double C[m][n]) {
	int i, j, p;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			C[i][j] = 0.0;
			for (p = 0; p < k; ++p) {
				C[i][j] += A[i][p]*B[p][j];
			}
		}
	}
}

// C = A^T * B
void MatTransposeMatMult_(const int m, const int k, const int n, \
	double A[k][m], double B[k][n], double C[m][n]) {
	int i, j, p;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			C[i][j] = 0.0;
			for (p = 0; p < k; ++p) {
				C[i][j] += A[p][i]*B[p][j];
			}
		}
	}
}

// C = A * B^T
void MatMatTransposeMult_(const int m, const int k, const int n, \
	double A[m][k], double B[n][k], double C[m][n]) {
	int i, j, p;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++ j) {
			C[i][j] = 0.0;
			for (p = 0; p < k; ++p) {
				C[i][j] += A[i][p]*B[j][p];
			}
		}
	}
}

// C = A + B
void MatMatAdd_(const int m, const int n, double A[m][n], double B[m][n], double C[m][n]) {
	int i, j;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

