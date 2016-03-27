#include "face.h"
#include <stdlib.h>


int FaceNum(const int nsd, const int nen) {
    int n;
    if (nsd == 2) {
        if (nen == 3) n = 3;
        else if (nen == 4) n = 4;
    } else if (nsd == 3) {
        if (nen == 4) n = 4;
        else if (nen == 8) n = 6;
    }
    return n;
}

int FaceNumNodes(const int nsd, const int nen) {
    int n;
    if (nsd == 2) {
        if (nen == 3 || nen == 4) n = 2;
    }
    else if (nsd == 3) {
        if (nen == 4) n = 3;
        else if (nen == 8) n = 4;
    }
    return n;
}


void FaceNodes(const int num_nodes, double p[num_nodes], const int face, \
    const int nsd, const int nen) {
    int i3[] = {2, 3, 1};
    int i4[] = {2, 3, 4, 1};
    if (nsd == 2) {
        if (nen == 3) {
            p[0] = face;
            p[1] = i3[face];
        } else if (nen == 4) {
            p[0] = face;
            p[1] = i4[face];
        }
    } else if (nsd == 3) {
        if (nen == 4) {
            if (face == 1) {
                p[0] = 1;
                p[1] = 2;
                p[2] = 3;
            } else if (face == 2) {
                p[0] = 1;
                p[1] = 4;
                p[2] = 2;
            } else if (face == 3) {
                p[0] = 2;
                p[1] = 4;
                p[2] = 3;
            } else if (face == 4) {
                p[0] = 3;
                p[1] = 4;
                p[2] = 1;
            }
        } else if (nen == 8) {
            if (face == 1) {
                p[0] = 1;
                p[1] = 2;
                p[2] = 3;
                p[3] = 4;
            } else if (face == 2) {
                p[0] = 5;
                p[1] = 8;
                p[2] = 7;
                p[3] = 6;
            } else if (face == 3) {
                p[0] = 1;
                p[1] = 5;
                p[2] = 6;
                p[3] = 2;
            } else if (face == 4) {
                p[0] = 2;
                p[1] = 6;
                p[2] = 7;
                p[3] = 3;
            } else if (face == 5) {
                p[0] = 3;
                p[1] = 7;
                p[2] = 8;
                p[3] = 4;
            } else if (face == 6) {
                p[0] = 4;
                p[1] = 8;
                p[2] = 5;
                p[3] = 1;
            }
        }
    }
}
