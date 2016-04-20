#ifndef EXTERNAL_FORCE_H
#define EXTERNAL_FORCE_H

int ExternalPressure(Vec* dofs, Vec* fglo);
int ExternalTraction(Vec* fglo);
int ExternalGravity(Vec* fglo);

#endif