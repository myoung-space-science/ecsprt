#ifndef RHS_H
#define RHS_H

#include <petsc.h>

extern const char *RHSTypes[];

typedef enum {
  RHS_CONSTANT,
  RHS_SINUSOIDAL,
  RHS_FULL,
} RHSType;

extern PetscErrorCode ComputeConstantRHS(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeSinusoidalRHS2D(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeSinusoidalRHS3D(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeFullRHS2D(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeFullRHS3D(KSP ksp, Vec b, void *opts);

typedef PetscErrorCode (*RHSFunc)(KSP ksp, Vec b, void *opts);

#endif // RHS_H
