#ifndef EXTERNAL_FORCE_H
#define EXTERNAL_FORCE_H

void ExternalPressure(int nsd, int nn, int nel, int nen, double dofs[nn*nsd+nel], double coords[nn][nsd], int connect[nel][nen], \
	int load_size, int load_num[load_size][2], int load_type, double load_val[load_size][load_type], double fglo[nn*nsd+nel]);
void ExternalTraction(int nsd, int nn, int nel, int nen, double coords[nn][nsd], int connect[nel][nen], \
	int load_size, int load_num[load_size][2], int load_type, double load_val[load_size][load_type], double fglo[nn*nsd+nel]);
void ExternalGravity(int nsd, int nn, int nel, int nen, double coords[nn][nsd], int connect[nel][nen], \
	double rho, double gravity[nsd], double fglo[nn*nsd+nel]);

#endif