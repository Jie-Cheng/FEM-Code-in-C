#include "helpers.h"
#include <stdlib.h>

extern void DGETRF_(int* M, int *N, double* A, int* LDA, int* IPIV, int* INFO);
extern void DGETRI_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
extern void DGEMM_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, \
	double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);

// Find the inverse of a matrix.
// A is a vectorized matrix, whose original size is N by N.
void Inverse(double* A, int N) {
	int IPIV[N];
    int LWORK = 2*N;
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

// Determinant of 2X2 or 3X3 matrix
double Determinant(const double** a, const int n) {
	double det;
	if (n == 2) {
		det = a[0][0]*a[1][1] - a[0][1]*a[1][0];
	} else if (n == 3) {
		det = a[0][0]*a[1][1]*a[2][2] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] \
		+  a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0];
	}
	return det;
}

// Calculates C = alpha*A*B + beta*C
// A[m][k] B[k][n] C[m][n] are originally 2d row-major matrices
void MatMul(int m, int n, int k, double alpha, double beta, double* A, double* B, double* C){
	// Change to column-major
	double *AT, *BT, *CT;
	int i, j, pos1, pos2;
	int lda = 2*m;
	int ldb = 2*k;
	int ldc = 2*m;
	AT = calloc(lda*k, sizeof(double));
	BT = calloc(ldb*n, sizeof(double));
	CT = calloc(ldc*n, sizeof(double));
	for (i = 0; i < m; ++i) {
		for (j = 0; j < k; ++j) {
			pos1 = k*i + j;
			pos2 = lda*j + i;
			AT[pos2] = A[pos1];
		}
	}
	for (i = 0; i < k; ++i) {
		for (j = 0; j < n; ++j) {
			pos1 = n*i + j;
			pos2 = ldb*j + i;
			BT[pos2] = B[pos1];
		}
	}
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			pos1 = n*i + j;
			pos2 = ldc*j + i;
			CT[pos2] = CT[pos1];
		}
	}
	DGEMM_("N", "N", &m, &n, &k, &alpha, AT, &lda, BT, &ldb, &beta, CT, &ldc);
	// Change the result to row-major
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			pos1 = n*i + j;
			pos2 = ldc*j + i;
			C[pos1] = CT[pos2];
		}
	}
	free(AT);
	free(BT);
	free(CT);
}