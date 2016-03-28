#ifndef HELPERS_H
#define HELPERS_H

void Inverse(double* A, int N);
void Cross3d(double* a, double* b, double* c);
double Determinant(const double** a, const int n);
void MatMul(int m, int n, int k, double alpha, double beta, double* A, double* B, double* C);

#endif