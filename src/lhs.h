#ifndef LHS_H
#define LHS_H

#include <petsc.h>
#include "constants.h"

extern const char *LHSTypes[];

typedef enum {
  LHS_IDENTITY,
  LHS_LAPLACIAN,
  LHS_FULL,
} LHSType;

extern PetscErrorCode ComputeIdentityLHS(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeLaplacianLHS(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeFullLHS(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeLHSEigenvalues(KSP ksp, void *opts);

typedef PetscErrorCode (*LHSFunc)(KSP ksp, Mat J, Mat A, void *opts);

#endif // LHS_H
