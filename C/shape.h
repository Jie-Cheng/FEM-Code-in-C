#ifndef SHAPE_H
#define SHAPE_H

void ShapeFun(const int nsd, double xi[nsd], const int nen, double p[nen]);
void ShapeDer(const int nsd, double xi[nsd], const int nen, double p[nen][nsd]);

#endif
