#ifndef MATERIAL_H
#define MATERIAL_H

void MaterialStiffness(int nsd, double intcoord[nsd], double F[nsd][nsd], \
	double pressure, int materialtype, double materialprops[5], double mstiff[nsd][nsd][nsd][nsd]);
void KirchhoffStress(int nsd, double intcoord[nsd], double F[nsd][nsd], double pressure, \
	int materialtype, double materialprops[5], double stress[nsd][nsd]);

#endif