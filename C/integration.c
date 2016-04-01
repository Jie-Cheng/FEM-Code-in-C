#include "integration.h"
#include <math.h>


int IntNum(const int nsd, const int nen, const int isreduced) {
    // The number of integration points in an element
    int n;
    if (!isreduced) {
        if (nsd == 1) {
            if (nen == 2) {
                n = 2;
            }
        } else if (nsd == 2) {
            if (nen == 3) {
                n = 1;
            } else if (nen == 4) {
                n = 4;
            }
        } else if (nsd == 3) {
            if (nen == 4) {
                n = 1;
            } else if (nen == 8) {
                n = 8;
            }
        }
    } else {
        if (nsd == 1) {
            if (nen == 2) {
                n = 2;
            }
        } else if (nsd == 2) {
            if (nen == 3) {
                n = 1;
            } else if (nen == 3) {
                n = 1;
            }
        } else if (nsd == 3) {
            if (nen == 4) {
                n = 1;
            } else if (nen == 8) {
                n = 1;
            }
        }
    }
    return n;
}


void IntWeights(const int nsd, const int nen, const int npt, double p[npt]) {
    // The weights of integration points
    if (nsd == 1) {
        if (npt == 2) {
            p[0] = 1.0;
            p[1] = 1.0;
        }
    } else if (nsd == 2) {
        if (nen == 3) {
            if (npt == 1) {
                p[0] = 0.5;
            }
        } else if (nen == 4) {
            if (npt == 1) {
                p[0] = 4.0;
            } else if (npt == 4) {
                p[0] = 1.0;
                p[1] = 1.0;
                p[2] = 1.0;
                p[3] = 1.0;
            }
        }
    } else if (nsd == 3) {
        if (nen == 4) {
            if (npt == 1) {
                p[0] = 1/6.0;
            }
        } else if (nen == 8) {
            if (npt == 1) {
                p[0] = 8.0;
            } else if (npt == 8) {
                p[0] = 1.0;
                p[1] = 1.0;
                p[2] = 1.0;
                p[3] = 1.0;
                p[4] = 1.0;
                p[5] = 1.0;
                p[6] = 1.0;
                p[7] = 1.0;
            }
        }
    }
}


void IntPoints(const int nsd, const int nen, const int npt, double p[npt][nsd]) {
    // The coordinates of the integration points
    int i, j;
    for (i = 0; i < npt; ++i) {
        for (j = 0; j < nsd; ++j) {
            p[i][j] = 0.0;
        }
    }
    double temp = 1/sqrt(3.0);
    if (nsd == 1) {
        if (npt == 2) {
            p[0][0] = -temp;
            p[1][0] = temp;
        }
    } else if (nsd == 2) {
        if (nen == 3) {
            if (npt == 1) {
                p[0][0] = 1/3.0;
                p[0][1] = 1/3.0;
            }
        } else if (nen == 4) {
            if (npt == 1) {
                p[0][0] = 0.0;
                p[0][1] = 0.0;
            } else if (npt == 4) {
                p[0][0] = -temp;
                p[0][1] = -temp;
                p[1][0] = temp;
                p[1][1] = -temp;
                p[2][0] = -temp;
                p[2][1] = temp;
                p[3][0] = temp;
                p[3][1] = temp;
            }
        }
    } else if (nsd == 3) {
        if (nen == 4) {
            if (npt == 1) {
                p[0][0] = 0.25;
                p[0][1] = 0.25;
                p[0][2] = 0.25;
            }
        } else if (nen == 8) {
            if (npt == 1) {
                p[0][0] = 0.0;
                p[0][1] = 0.0;
                p[0][2] = 0.0;
            } else if (npt == 8) {
                p[0][0] = -1.0;
                p[0][1] = -1.0;
                p[0][2] = -1.0;
                p[1][0] = 1.0;
                p[1][1] = -1.0;
                p[1][2] = -1.0;
                p[2][0] = 1.0;
                p[2][1] = 1.0;
                p[2][2] = -1.0;
                p[3][0] = -1.0;
                p[3][1] = 1.0;
                p[3][2] = -1.0;
                p[4][0] = -1.0;
                p[4][1] = -1.0;
                p[4][2] = 1.0;
                p[5][0] = 1.0;
                p[5][1] = -1.0;
                p[5][2] = 1.0;
                p[6][0] = 1.0;
                p[6][1] = 1.0;
                p[6][2] = 1.0;
                p[7][0] = -1.0;
                p[7][1] = 1.0;
                p[7][2] = 1.0;
            }
        }
    }
}
