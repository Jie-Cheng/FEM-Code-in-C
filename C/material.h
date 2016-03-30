#ifndef MATERIAL_H
#define MATERIAL_H

void MaterialStiffness(int nsd, double intcoord[nsd], double F[nsd][nsd], \
	double pressure, int materialtype, double* materialprops, double mstiff[nsd][nsd][nsd][nsd]);
void KirchhoffStress(int nsd, double intcoord[nsd], double F[nsd][nsd], double pressure, \
	int materialtype, double* materialprops, double stress[nsd][nsd]);

#endif