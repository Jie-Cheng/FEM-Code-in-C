#include "helpers.h"
#include <math.h>

extern void DGETRF_(int* M, int *N, double* A, int* LDA, int* IPIV, int* INFO);
extern void DGETRI_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
extern void DGEMM_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, \
	double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);

// Find the inverse of a matrix.
// A is a vectorized matrix, whose original size is N by N.
void Inverse(int N, double* A) {
	int IPIV[N];
    int LWORK = 2*N;
    double WORK[LWORK];
    int INFO;

    DGETRF_(&N,&N,A,&N,IPIV,&INFO);
    DGETRI_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
}


// Find the cross product of two vectors.
void Cross3d(double a[3], double b[3], double c[3]) {
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

// Determinant of 2X2 or 3X3 matrix
double Determinant(const int n, double a[n][n]) {
	// Assuming a[n][n]
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
void MatMul(char transa, char transb, int m, int n, int k, double alpha, double beta, double* A, double* B, double* C){
	int lda, ldb, ldc;
	int i, j, pos1, pos2;
	int ka, kb;

	ldc = 2*m;

	if (transa == 'N') {
		lda = 2*m;
		ka = k;
	} else if (transa == 'T') {
		lda = 2*k;
		ka = m;
	}

	if (transb == 'N') {
		ldb = 2*k;
		kb = n;
	} else if (transb == 'T') {
		ldb = 2*n;
		kb = k;
	}

	double OA[lda*ka], OB[ldb*kb], OC[ldc*n];

	if (transa == 'N') {
		// A is m by k
		// Change order because of C/FORTRAN difference
		for (i = 0; i < m; ++i) {
			for (j = 0; j < k; ++j) {
				pos1 = k*i + j;
				pos2 = lda*j + i;
				OA[pos2] = A[pos1];
			}
		}
	} else if (transa == 'T') {
		// A is k by m
		// Do not change
		for (i = 0; i < k; ++i) {
			for (j = 0; j < m; ++j) {
				pos1 = m*i + j;
				pos2 = lda*j + i;
				OA[pos2] = A[pos1];
			}
		}
	}

	if (transb == 'N') {
		// B is k by n
		for (i = 0; i < k; ++i) {
			for (j = 0; j < n; ++j) {
				pos1 = n*i + j;
				pos2 = ldb*j + i;
				OB[pos2] = B[pos1];
			}
		}
	} else if (transb == 'T') {
		// B is n by k
		for (i = 0; i < n; ++i) {
			for (j = 0; j < k; ++j) {
				pos1 = k*i + j;
				pos2 = ldb*j + i;
				OB[pos2] = B[pos1];
			}
		}
	}

	if (fabs(beta) > 1e-5) {
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				pos1 = n*i + j;
				pos2 = ldc*j + i;
				OC[pos2] = C[pos1];
			}
		}
	}
	

	DGEMM_(&transa, &transb, &m, &n, &k, &alpha, OA, &lda, OB, &ldb, &beta, OC, &ldc);

	// Change the result to row-major
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			pos1 = n*i + j;
			pos2 = ldc*j + i;
			C[pos1] = OC[pos2];
		}
	}
}

// Dot product
double DotProduct(const int n, double* a, double* b) {
	// Assuming a[n], b[n]
	int i;
	double result = 0.0;
	for (i = 0; i < n; ++i) {
		result += a[i]*b[i];
	}
	return result;
}