#ifndef VELOCITIES_H
#define VELOCITIES_H

#include <petsc.h>
#include "ecsprt.h"

extern PetscErrorCode NormalVelocities(Context *ctx);
extern PetscErrorCode BorisMoverBB(PetscReal dt, Context *ctx);
extern PetscErrorCode BorisMoverBz(PetscReal dt, Context *ctx);

#endif // VELOCITIES_H