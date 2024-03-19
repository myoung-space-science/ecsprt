#ifndef VELOCITIES_H
#define VELOCITIES_H

#include <petsc.h>
#include "ecsprt.h"
#include "pic.h"

extern PetscErrorCode NormalVelocities(Context *ctx);
extern PetscErrorCode InitializeVelocities(VDistType VDistType, Context *ctx);
extern PetscErrorCode BorisMoverBB(PetscReal dt, Context *ctx);
extern PetscErrorCode BorisMoverBz(PetscReal dt, Context *ctx);
extern PetscErrorCode UpdateVelocities(PetscReal dt, Context *ctx);

#endif // VELOCITIES_H