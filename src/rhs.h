#ifndef RHS_H
#define RHS_H

#include <petsc.h>

extern PetscErrorCode ComputeConstantRHS(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeSinusoidalRHS2D(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeSinusoidalRHS3D(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeFullRHS2D(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeFullRHS3D(KSP ksp, Vec b, void *opts);

#endif // RHS_H
