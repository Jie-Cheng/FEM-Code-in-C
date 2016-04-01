#ifndef HELPERS_H
#define HELPERS_H
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))

void Inverse(int N, double* A);
void Cross3d(double a[3], double b[3], double c[3]);
double Determinant(const int n, double a[n][n]);
void MatMul(char transa, char transb, int m, int n, int k, double alpha, double beta, double* A, double* B, double* C);
double DotProduct(const int n, double* a, double* b);

#endif