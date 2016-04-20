#ifndef MATERIAL_H
#define MATERIAL_H

void KirchhoffStress(double intcoord[nsd], double F[nsd][nsd], double pressure, double stress[nsd][nsd]);
void MaterialStiffness(double intcoord[nsd], double F[nsd][nsd], double pressure, double mstiff[nsd][nsd][nsd][nsd]);

#endif