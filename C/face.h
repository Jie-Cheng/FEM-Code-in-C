#ifndef FACE_H
#define FACE_H

int FaceNum(const int nsd, const int nen);
int FaceNumNodes(const int nsd, const int nen);
void FaceNodes(const int nsd, const int nen, const int face, \
    const int num_nodes, double p[num_nodes]);

#endif
