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
extern PetscErrorCode ComputeSinusoidalRHS(KSP ksp, Vec b, void *opts);
extern PetscErrorCode ComputeFullRHS(KSP ksp, Vec b, void *opts);

typedef PetscErrorCode (*RHSFunc)(KSP ksp, Vec b, void *opts);

#endif // RHS_H
