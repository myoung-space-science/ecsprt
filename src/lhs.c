#include <petsc.h>
#include <slepceps.h>
#include "constants.h"
#include "ecsprt.h"
#include "lhs.h"

PetscErrorCode ComputeIdentityLHS(KSP ksp, Mat J, Mat A, void *opts)
{
  Context      *ctx=(Context *)opts;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(MatZeroEntries(A));
  PetscCall(MatShift(A, 1.0));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeLaplacianLHS(KSP ksp, Mat J, Mat A, void *opts)
{
  Context      *ctx=(Context *)opts;
  PetscReal     dx=ctx->grid.dx;
  PetscReal     dy=ctx->grid.dy;
  PetscReal     dz=ctx->grid.dz;
  DM            dm;
  PetscInt      i0, j0, k0;
  PetscInt      ni, nj, nk;
  PetscInt      i, j, k;
  PetscScalar   vijk=1.0;
  PetscScalar   vpjk=0.0, vmjk=0.0, vipk=0.0, vimk=0.0, vijp=0.0, vijm=0.0;
  MatStencil    row;
  PetscScalar   vals[ctx->potential.stencilSize];
  MatStencil    cols[ctx->potential.stencilSize];
  MatNullSpace  nullspace;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Assign the star-stencil coefficients.
  vpjk =  1.0 / (dx*dx);
  vmjk =  1.0 / (dx*dx);
  vipk =  1.0 / (dy*dy);
  vimk =  1.0 / (dy*dy);
  vijp =  1.0 / (dz*dz);
  vijm =  1.0 / (dz*dz);

  // Assign the diagonal coefficient.
  vijk = -(vpjk + vipk + vijp + vmjk + vimk + vijm);

  // Get the DM associated with the KSP.
  PetscCall(KSPGetDM(ksp, &dm));

  // Get this processor's indices.
  PetscCall(DMDAGetCorners(dm, &i0, &j0, &k0, &ni, &nj, &nk));

  // Loop over grid points. [DEV] Assume periodic BC.
  for (k=k0; k<k0+nk; k++) {
    for (j=j0; j<j0+nj; j++) {
      for (i=i0; i<i0+ni; i++) {
        row.i = i; row.j = j; row.k = k;
        // Assign the value at node (i+1, j, k)
        vals[0] = vpjk;
        cols[0].i = i+1;
        cols[0].j = j;
        cols[0].k = k;
        // Assign the value at node (i-1, j, k)
        vals[1] = vmjk;
        cols[1].i = i-1;
        cols[1].j = j;
        cols[1].k = k;
        // Assign the value at node (i, j+1, k)
        vals[2] = vipk;
        cols[2].i = i;
        cols[2].j = j+1;
        cols[2].k = k;
        // Assign the value at node (i, j-1, k)
        vals[3] = vimk;
        cols[3].i = i;
        cols[3].j = j-1;
        cols[3].k = k;
        // Assign the value at node (i, j, k+1)
        vals[4] = vijp;
        cols[4].i = i;
        cols[4].j = j;
        cols[4].k = k+1;
        // Assign the value at node (i, j, k-1)
        vals[5] = vijm;
        cols[5].i = i;
        cols[5].j = j;
        cols[5].k = k-1;
        // Assign the value at node (i, j, k)
        vals[6] = vijk;
        cols[6].i = i;
        cols[6].j = j;
        cols[6].k = k;
        PetscCall(MatSetValuesStencil(A, 1, &row, ctx->potential.stencilSize, cols, vals, INSERT_VALUES));
      }
    }
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
  PetscCall(MatSetNullSpace(A, nullspace));
  PetscCall(MatNullSpaceDestroy(&nullspace));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeFullLHS(KSP ksp, Mat J, Mat A, void *opts)
{
  Context        *ctx=(Context *)opts;
  DM              fluidDM=ctx->fluidDM;
  PetscReal       kappa=ctx->electrons.kappa;
  PetscReal       dx=ctx->grid.dx;
  PetscReal       dy=ctx->grid.dy;
  PetscReal       dz=ctx->grid.dz;
  PetscReal       hxx, hyy, hzz, hod;
  PetscReal       scale=ctx->potential.scale;
  PetscReal       sxx, syy, szz, sod;
  Vec             moments;
  FluidNode    ***fluid;
  DM              dm;
  PetscInt        i0, j0, k0;
  PetscInt        ni, nj, nk;
  PetscInt        i, j, k;
  PetscReal       nijk, npjk, nmjk, nipk, nimk, nijp, nijm;
  Vec             F;
  PetscReal    ***f;
  MatStencil      row;
  PetscReal       vals[ctx->potential.stencilSize];
  MatStencil      cols[ctx->potential.stencilSize];
  MatNullSpace    nullspace;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Compute geometric scale factors for stencil values.
  hxx = 1.0 / (2.0 * dx*dx);
  hyy = 1.0 / (2.0 * dy*dy);
  hzz = 1.0 / (2.0 * dz*dz);
  hod = 1.0 / (8.0 * dy*dx);

  // Compute coefficient scale factors;
  sxx = scale * hxx;
  syy = scale * hyy;
  szz = scale * (1 + kappa*kappa)*hzz;
  sod = scale * kappa*hod;

  // Get density and flux arrays.
  PetscCall(DMGetLocalVector(fluidDM, &moments));
  PetscCall(DMGlobalToLocalBegin(fluidDM, ctx->moments, INSERT_VALUES, moments));
  PetscCall(DMGlobalToLocalEnd(fluidDM, ctx->moments, INSERT_VALUES, moments));
  PetscCall(DMDAVecGetArray(fluidDM, moments, &fluid));

  // Get the DM associated with the KSP.
  PetscCall(KSPGetDM(ksp, &dm));

  // Get a local fluid of zeros that shares data with the KSP DM.
  PetscCall(DMGetLocalVector(dm, &F));
  PetscCall(VecZeroEntries(F));
  PetscCall(DMDAVecGetArray(dm, F, &f));

  // Get this processor's indices.
  PetscCall(DMDAGetCorners(dm, &i0, &j0, &k0, &ni, &nj, &nk));

  // Loop over grid points.
  for (k=k0; k<k0+nk; k++) {
    for (j=j0; j<j0+nj; j++) {
      for (i=i0; i<i0+ni; i++) {

        // Assign density values.
        nijk = fluid[k][j][i].n;
        nmjk = fluid[k][j][i-1].n;
        npjk = fluid[k][j][i+1].n;
        nimk = fluid[k][j-1][i].n;
        nipk = fluid[k][j+1][i].n;
        nijm = fluid[k-1][j][i].n;
        nijp = fluid[k+1][j][i].n;

        // Assign the x-y corner coefficients
        f[k][j+1][i+1] = sod*(nipk - npjk);
        f[k][j-1][i+1] = sod*(npjk - nimk);
        f[k][j+1][i-1] = sod*(nmjk - nipk);
        f[k][j-1][i-1] = sod*(nimk - nmjk);

        // Assign the star-stencil coefficients
        f[k][j][i+1] =  sxx*(npjk + nijk) + sod*(nipk - nimk);
        f[k][j][i-1] =  sxx*(nijk + nmjk) - sod*(nipk - nimk);
        f[k][j+1][i] =  syy*(nipk + nijk) - sod*(npjk - nmjk);
        f[k][j-1][i] =  syy*(nijk + nimk) + sod*(npjk - nmjk);
        f[k+1][j][i] =  szz*(nijp + nijk);
        f[k-1][j][i] =  szz*(nijk + nijm);

        // Assign the diagonal coefficient
        f[k][j][i] = -(sxx*(npjk + 2*nijk + nmjk) + syy*(nipk + 2*nijk + nimk) + szz*(nijp + 2*nijk + nijm));

        // Assign the value at node (i+1, j, k).
        vals[0] = f[k][j][i+1];
        cols[0].i = i+1;
        cols[0].j = j;
        cols[0].k = k;
        // Assign the value at node (i-1, j, k).
        vals[1] = f[k][j][i-1];
        cols[1].i = i-1;
        cols[1].j = j;
        cols[1].k = k;
        // Assign the value at node (i, j+1, k).
        vals[2] = f[k][j+1][i];
        cols[2].i = i;
        cols[2].j = j+1;
        cols[2].k = k;
        // Assign the value at node (i, j-1, k).
        vals[3] = f[k][j-1][i];
        cols[3].i = i;
        cols[3].j = j-1;
        cols[3].k = k;
        // Assign the value at node (i, j, k+1).
        vals[4] = f[k+1][j][i];
        cols[4].i = i;
        cols[4].j = j;
        cols[4].k = k+1;
        // Assign the value at node (i, j, k-1).
        vals[5] = f[k-1][j][i];
        cols[5].i = i;
        cols[5].j = j;
        cols[5].k = k-1;
        // Assign the value at node (i+1, j+1, k).
        vals[6] = f[k][j+1][i+1];
        cols[6].i = i+1;
        cols[6].j = j+1;
        cols[6].k = k;
        // Assign the value at node (i+1, j-1, k).
        vals[7] = f[k][j-1][i+1];
        cols[7].i = i+1;
        cols[7].j = j-1;
        cols[7].k = k;
        // Assign the value at node (i-1, j+1, k).
        vals[8] = f[k][j+1][i-1];
        cols[8].i = i-1;
        cols[8].j = j+1;
        cols[8].k = k;
        // Assign the value at node (i-1, j-1, k).
        vals[9] = f[k][j-1][i-1];
        cols[9].i = i-1;
        cols[9].j = j-1;
        cols[9].k = k;
        // Assign the value at node (i, j, k).
        vals[10] = f[k][j][i];
        cols[10].i = i;
        cols[10].j = j;
        cols[10].k = k;
        row.i = i; row.j = j; row.k = k;
        PetscCall(MatSetValuesStencil(A, 1, &row, ctx->potential.stencilSize, cols, vals, INSERT_VALUES));
      }
    }
  }

  // Restore the coefficient array and corresponding vector.
  PetscCall(DMDAVecRestoreArray(dm, F, &f));
  PetscCall(DMRestoreLocalVector(dm, &F));

  // Restore density array and corresponding vector.
  PetscCall(DMDAVecRestoreArray(fluidDM, moments, &fluid));
  PetscCall(DMRestoreLocalVector(fluidDM, &moments));

  // Assemble the distributed operator matrix.
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  // Remove the operator null space.
  PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
  PetscCall(MatSetNullSpace(A, nullspace));
  PetscCall(MatNullSpaceDestroy(&nullspace));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the eigenvalues of the operator matrix.

This is essentially a distilled version of ${SLEPC_DIR}/src/eps/tutorials/ex1.c
*/
PetscErrorCode ComputeLHSEigenvalues(KSP ksp, void *opts)
{
  Context     *ctx=(Context *)opts;
  EPS          eps;
  EPSType      type;
  Mat          A;
  PetscInt     i, its, maxit, nev, nconv;
  PetscReal    tol, error, re, im;
  PetscScalar  kr, ki;
  Vec          xr, xi;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  /* Set up eigenvalue solver. */
  PetscCall(EPSCreate(PETSC_COMM_WORLD, &eps));

  /* Get the LHS operator from the Krylov context. */
  PetscCall(KSPGetOperators(ksp, &A, NULL));

  /* Create vectors for holding the real and imaginary parts. */
  PetscCall(MatCreateVecs(A, NULL, &xr));
  PetscCall(MatCreateVecs(A, NULL, &xi));

  /* Associate the LHS operator with the eigenvalue solver. */
  PetscCall(EPSSetOperators(eps, A, NULL));

  /* Set up the problem type and read CLI options. */
  PetscCall(EPSSetProblemType(eps, EPS_HEP));
  PetscCall(EPSSetFromOptions(eps));

  /* Compute eigenvalues. */
  PetscCall(EPSSolve(eps));

  /* Print results. */
  PetscCall(EPSGetIterationNumber(eps, &its));
  ctx->log.status(" Number of iterations of the method: %d\n", its);
  PetscCall(EPSGetType(eps, &type));
  ctx->log.status(" Solution method: %s\n\n", type);
  PetscCall(EPSGetDimensions(eps, &nev, NULL, NULL));
  ctx->log.status(" Number of requested eigenvalues: %d\n", nev);
  PetscCall(EPSGetTolerances(eps, &tol, &maxit));
  ctx->log.status(" Stopping condition: tol=%.4g, maxit=%d\n", (double)tol, maxit);
  PetscCall(EPSGetConverged(eps, &nconv));
  ctx->log.status(" Number of converged eigenpairs: %d\n\n", nconv);
  if (nconv > 0) {
    ctx->log.status(
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");
    for (i=0; i<nconv; i++) {
      PetscCall(EPSGetEigenpair(eps, i, &kr, &ki, xr, xi));
      PetscCall(EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error));
#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) ctx->log.status(" %9f%+9fi %12g\n", (double)re, (double)im, (double)error);
      else ctx->log.status("   %12f       %12g\n", (double)re, (double)error);
    }
    ctx->log.status("\n");
  }

  /* Free memory. */
  PetscCall(EPSDestroy(&eps));
  PetscCall(VecDestroy(&xr));
  PetscCall(VecDestroy(&xi));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}