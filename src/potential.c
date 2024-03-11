#include <petsc.h>
#include "ecsprt.h"


/* Compute the electrostatic potential and store the result in the context. */
PetscErrorCode ComputePotential(KSP ksp, Context *ctx)
{
  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(KSPSolve(ksp, NULL, NULL));
  PetscCall(KSPGetSolution(ksp, &ctx->potential.solution));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Set the initial guess for the electrostatic potential. */
PetscErrorCode ComputeInitialPhi(KSP ksp, Vec phi, void *opts)
{
  // Note that `KSPSetComputeInitialGuess` requires this function signature.

  PetscFunctionBeginUser;

  PetscCall(VecSet(phi, 0.0));

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the forcing vector for the electrostatic-potential equation. */
PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *opts)
{
  Context      *ctx=(Context *)opts;

  PetscFunctionBeginUser;

  PetscCall(ctx->potential.rhs(ksp, b, opts));

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the operator matrix for the electrostatic-potential equation. */
PetscErrorCode ComputeLHS(KSP ksp, Mat J, Mat A, void *opts)
{
  Context      *ctx=(Context *)opts;

  PetscFunctionBeginUser;

  PetscCall(ctx->potential.lhs(ksp, J, A, opts));

  PetscFunctionReturn(PETSC_SUCCESS);
}


