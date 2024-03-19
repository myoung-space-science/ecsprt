#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "ecsprt.h"

extern PetscErrorCode Apply2DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts);
extern PetscErrorCode Apply3DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts);
extern PetscErrorCode ApplyBCAndMigrate(Context *ctx);

#endif // BOUNDARIES_H
