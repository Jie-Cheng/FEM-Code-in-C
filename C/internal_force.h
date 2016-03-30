#ifndef INTERNAL_FORCE_H
#define INTERNAL_FORCE_H

void InternalForce(int nsd, int nn, int nel, int nen, double dofs[nn*nsd+nel], \
	double coords[nn][nsd], int connect[nel][nen], int materialtype, \
	double* materialprops, double fglo[nn*nsd+nel]); 

#endif