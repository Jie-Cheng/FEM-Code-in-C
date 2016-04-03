#ifndef READ_H
#define READ_H

#ifdef __cplusplus
extern "C" {
#endif
    void ReadInput(int* mode, int* maxit, int* nsteps, int* nprint, int* isbinary, \
    double* inistep, double* adjust, double* tol, double* dt, double* damp, double* penalty, \
    int* materialtype, double materialprops[5], double gravity[3]);
    
    void ReadMesh(int* nsd, int* nn, int* nel, int* nen, double coords[*nn][*nsd], int connect[*nel][*nen], \
    int* bc_size, int bc_num[*bc_size][2], double bc_val[*bc_size], int* load_size, int* load_type, \
    int load_num[*load_size][2], double load_val[*load_size][*load_type], int share[*nn]);

#ifdef __cplusplus
}
#endif

#endif
