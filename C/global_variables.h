#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

extern int mode, maxit, nsteps, nprint, step, isbinary;
extern double inistep, adjust, tol, dt, damp, penalty;
extern int materialtype;
extern double materialprops[5];
extern double gravity[3];
extern int nsd, nen, nn, nel;
extern int** connect;
extern double** coords;
extern int** bc_num;
extern double* bc_val;
extern int** load_num;
extern double** load_val;
extern int* share;


#endif
