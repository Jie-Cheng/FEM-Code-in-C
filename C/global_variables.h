#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

int mode, maxit, nsteps, nprint, step, isbinary;
double inistep, adjust, tol, dt, damp, penalty;
int materialtype;
double materialprops[5];
double gravity[3];
int nsd, nen, nn, nel;
int** connect;
double** coords;
int** bc_num;
double* bc_val;
int** load_num;
double** load_val;
int* share;


#endif