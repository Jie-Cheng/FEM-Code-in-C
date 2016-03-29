#include <stdio.h>
#include "global_variables.h"

extern void ReadInput();
extern void ReadMesh();

int main() {
    ReadInput();
    ReadMesh();
    return 0;

}
