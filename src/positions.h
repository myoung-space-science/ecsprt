#ifndef POSITIONS_H
#define POSITIONS_H

#include <petsc.h>
#include "ecsprt.h"
#include "pic.h"

// Type to be used for particle distribution functions in Rejection.
typedef PetscErrorCode (*DistributionFunction)(PetscInt ndim, PetscReal r[], PetscReal *v, Context *ctx);

// Top-level distribution functions.
extern PetscErrorCode UniformCoordinates(PetscInt ndim, Context *ctx);
extern PetscErrorCode SobolDistribution(PetscInt ndim, Context *ctx);
extern PetscErrorCode NormalDistribution(PetscInt ndim, Context *ctx);
extern PetscErrorCode Rejection(PetscInt ndim, DistributionFunction density, Context *ctx);
extern PetscErrorCode InitializePositions(PetscInt ndim, PDistType PDistType, Context *ctx);
extern PetscErrorCode UpdatePositions(PetscInt ndim, PetscReal dt, Context *ctx);

// DistributionFunction implementations.
extern PetscErrorCode SinusoidalDistribution(PetscInt ndim, PetscReal r[], PetscReal *v, Context *ctx);

#endif // POSITIONS_H