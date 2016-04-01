#ifndef FACE_H
#define FACE_H

int FaceNum(const int nsd, const int nen);
int FaceNumNodes(const int nsd, const int nen);
void FaceNodes(const int nsd, const int nen, const int nfacenodes, const int face, int p[nfacenodes]);

#endif
