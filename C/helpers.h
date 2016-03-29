#ifndef HELPERS_H
#define HELPERS_H
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))
#define DELTA(x, y) (((x) == (y)) ? (1.0) : (0.0))
#define M_PI 3.14159265358979323846

void Inverse(double* A, int N);
void Cross3d(double* a, double* b, double* c);
double Determinant(const double** a, const int n);
void MatMul(char transa, char transb, int m, int n, int k, double alpha, double beta, double* A, double* B, double* C);
double DotProduct(int n, double a[n], double b[n]);

#endif