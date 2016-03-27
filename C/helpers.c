#include "helpers.h"

extern void DGETRF_(int* M, int *N, double* A, int* LDA, int* IPIV, int* INFO);
extern void DGETRI_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);

// Find the inverse of a matrix.
// A is a vectorized matrix, whose original size is N by N.
void Inverse(double* A, int N) {
	int IPIV[N];
    int LWORK = N*N;
    double WORK[LWORK];
    int INFO;

    DGETRF_(&N,&N,A,&N,IPIV,&INFO);
    DGETRI_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
}


// Find the cross product of two vectors.
// Assuming a, b, c are of size 3.
void Cross3d(double* a, double* b, double* c) {
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}