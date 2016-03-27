#ifndef SHAPE_H
#define SHAPE_H

void ShapeFun(const int nen, double p[nen], const int nsd, double xi[nsd]);
void ShapeDer(const int nsd, const int nen, double p[nsd][nen], double xi[nsd]);

#endif
