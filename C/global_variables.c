#include "global_variables.h"

int mode, maxit, nsteps, nprint, step, isbinary;
double inistep, adjust, tol, dt, damp, penalty;
int materialtype;
double materialprops[5];
double gravity[3];
int nsd, nen, nn, nel;
int** connect = 0;
double** coords = 0;
int bc_size;
int** bc_num = 0;
double* bc_val = 0;
int load_size;
int load_type; // The number of columns of the load_val, could be 1, 2, 3
int** load_num = 0;
double** load_val = 0;
int* share = 0;
int* nnz = 0;