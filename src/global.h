#ifndef GLOBAL_H
#define GLOBAL_H

#include <petsc.h>

/* Defined types of LHS operators for the potential equation. */
typedef enum {
  LHS_IDENTITY,
  LHS_LAPLACIAN,
  LHS_FULL,
} LHSType;

/* Array that will hold names of LHS operators. */
extern const char *LHSTypes[];

/* Defined types of RHS vectors for the potential equation. */
typedef enum {
  RHS_CONSTANT,
  RHS_SINUSOIDAL,
  RHS_FULL,
} RHSType;

/* Array that will hold names of RHS operators. */
extern const char *RHSTypes[];

/* Function types. */
typedef PetscErrorCode (*LHSFunc)(KSP ksp, Mat J, Mat A, void *opts);
typedef PetscErrorCode (*RHSFunc)(KSP ksp, Vec b, void *opts);
typedef PetscErrorCode (*BCFunc)(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts);
typedef PetscErrorCode (*CollectFunc)(void *opts);
typedef PetscErrorCode (*CollideFunc)(PetscReal dt, void *opts);
typedef PetscErrorCode (*PushFunc)(PetscReal dt, void *opts);

#endif // GLOBAL_H
