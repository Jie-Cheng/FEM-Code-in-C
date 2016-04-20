#ifndef HELPERS_H
#define HELPERS_H
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))

extern double i3[3][3];
extern double i2[2][2];

void Cross3d(double a[3], double b[3], double c[3]);
double Determinant(const int n, double a[n][n]);
void InvertDirect(const int n, double A[n][n], double B[n][n]);
double DotProduct(const int n, double a[n], double b[n]);
void VecScale_(const double alpha, const int n, double v[n]);
void MatScale_(const double alpha, const int m, const int n, double A[m][n]);
void MatMult_(const int m, const int n, double A[m][n], double x[n], double b[m]);
void MatMultTranspose_(const int m, const int n, double A[m][n], double x[m], double b[n]);
void MatMatMult_(const int m, const int k, const int n, \
	double A[m][k], double B[k][n], double C[m][n]);
void MatTransposeMatMult_(const int m, const int k, const int n, \
	double A[k][m], double B[k][n], double C[m][n]);
void MatMatTransposeMult_(const int m, const int k, const int n, \
	double A[m][k], double B[n][k], double C[m][n]);
void MatMatAdd_(const int m, const int n, double A[m][n], double B[m][n], double C[m][n]);

#endif