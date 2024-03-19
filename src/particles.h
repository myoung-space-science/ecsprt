#ifndef PARTICLES_H
#define PARTICLES_H

#include <petsc.h>
#include "ecsprt.h"
#include "pic.h"

extern PetscErrorCode InitializePositions(PDistType PDistType, Context *ctx);
extern PetscErrorCode InitializeVelocities(VDistType VDistType, Context *ctx);
extern PetscErrorCode UpdateVelocities(PetscReal dt, Context *ctx);
extern PetscErrorCode UpdatePositions(PetscReal dt, Context *ctx);

#endif // PARTICLES_H