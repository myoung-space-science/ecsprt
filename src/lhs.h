#ifndef LHS_H
#define LHS_H

#include <petsc.h>
#include "constants.h"

extern PetscErrorCode ComputeIdentityLHS(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeLaplacianLHS2D(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeLaplacianLHS3D(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeFullLHS2D(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeFullLHS3D(KSP ksp, Mat J, Mat A, void *opts);
extern PetscErrorCode ComputeLHSEigenvalues(KSP ksp, void *opts);

#endif // LHS_H
