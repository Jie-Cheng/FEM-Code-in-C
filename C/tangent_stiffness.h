#ifndef TANGENT_STIFFNESS_H
#define TANGENT_STIFFNESS_H

void InternalTangent(int nsd, int nn, int nel, int nen, double dofs[nn*nsd+nel], \
	double coords[nn][nsd], int connect[nel][nen], int materialtype, \
	double materialprops[5], double kglo[nn*nsd+nel][nn*nsd+nel]);

#endif