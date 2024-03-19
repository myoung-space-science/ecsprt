#ifndef PARTICLES_H
#define PARTICLES_H

#include <petsc.h>
#include "ecsprt.h"
#include "pic.h"

extern PetscErrorCode Apply2DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts);
extern PetscErrorCode Apply3DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts);
extern PetscErrorCode InitializePositions(PDistType PDistType, Context *ctx);
extern PetscErrorCode InitializeVelocities(VDistType VDistType, Context *ctx);
extern PetscErrorCode BorisMover3D(PetscReal dt, Context *ctx);
extern PetscErrorCode BorisMoverBz(PetscReal dt, Context *ctx);
extern PetscErrorCode ComputeCollisions(PetscReal dt, Context *ctx);
extern PetscErrorCode UpdateVelocities(PetscReal dt, Context *ctx);
extern PetscErrorCode UpdatePositions(PetscReal dt, Context *ctx);

#endif // PARTICLES_H