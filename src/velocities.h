#ifndef VELOCITIES_H
#define VELOCITIES_H

#include <petsc.h>
#include "ecsprt.h"
#include "pic.h"

extern PetscErrorCode NormalVelocities(PetscInt ndim, Context *ctx);
extern PetscErrorCode InitializeVelocities(PetscInt ndim, VDistType vDistType, Context *ctx);
extern PetscErrorCode BorisMoverBB(PetscReal dt, Context *ctx);
extern PetscErrorCode BorisMoverBz2D(PetscReal dt, void *opts);
extern PetscErrorCode BorisMoverBz3D(PetscReal dt, void *opts);
extern PetscErrorCode UpdateVelocities(PetscInt ndim, PetscReal dt, Context *ctx);

#endif // VELOCITIES_H