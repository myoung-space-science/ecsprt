#ifndef PARTICLES_H
#define PARTICLES_H

#include <petsc.h>
#include "ecsprt.h"
#include "pic.h"

extern PetscErrorCode InitializePositions(DensityType densityType, Context *ctx);
extern PetscErrorCode InitializeVelocities(Context *ctx);
extern PetscErrorCode CollectFluidMoments(Context *ctx);
extern PetscErrorCode BorisMover3D(PetscReal dt, Context *ctx);
extern PetscErrorCode BorisMoverBz(PetscReal dt, Context *ctx);
extern PetscErrorCode ComputeCollisions(PetscReal dt, Context *ctx);
extern PetscErrorCode UpdateVelocities(PetscReal dt, Context *ctx);
extern PetscErrorCode UpdatePositions(PetscReal dt, Context *ctx);

#endif // PARTICLES_H