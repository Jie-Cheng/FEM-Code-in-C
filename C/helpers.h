#ifndef HELPERS_H
#define HELPERS_H
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))

void Inverse(double* A, int N);
void Cross3d(double* a, double* b, double* c);
double Determinant(const int n, const double a[n][n]);
void MatMul(char transa, char transb, int m, int n, int k, double alpha, double beta, double* A, double* B, double* C);
double DotProduct(const int n, const double a[n], const double b[n]);

#endif