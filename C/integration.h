#ifndef INTEGRATION_H
#define INTEGRATION_H

int IntNum(const int nsd, const int nen);
void IntWeights(const int nsd, const int nen, const int npt, double p[npt]);
void IntPoints(const int nsd, const int nen, const int npt, double p[npt][nsd]);

#endif
