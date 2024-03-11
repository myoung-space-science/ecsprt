#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <petsc.h>
#include "ecsprt.h"

extern PetscErrorCode ComputePotential(KSP ksp, Context *ctx);
extern PetscErrorCode ComputeInitialPhi(KSP ksp, Vec phi, void *opts);
extern PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeLHS(KSP ksp, Mat J, Mat A, void *opts);

#endif // POTENTIAL_H