#include <petsc.h>
#include "vector-math.h"
#include "hybrid.h"
#include "rhs.h"
#include "particles.h"

PetscErrorCode ComputeConstantRHS(KSP ksp, Vec b, void *user)
{
  Context   *ctx=(Context *)user;
  PetscReal  val;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Set the RHS vector equal to the background density.
  val = ctx->plasma.n0 * ctx->potential.scale;
  PetscCall(VecSet(b, val));

  // Store the RHS vector in the problem context.
  ctx->potential.forcing = b;

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeSinusoidalRHS(KSP ksp, Vec b, void *user)
{
  Context     *ctx=(Context *)user;
  PetscReal    scale=ctx->potential.scale;
  PetscReal    n0=ctx->plasma.n0;
  PetscReal    dx=ctx->grid.dx;
  PetscReal    dy=ctx->grid.dy;
  PetscReal    dz=ctx->grid.dz;
  DM           dm;
  PetscReal ***rhs;
  PetscInt     i0, j0, k0;
  PetscInt     ni, nj, nk;
  PetscInt     i, j, k;
  PetscReal    x, y, z;
  PetscReal    Cx, Cy, Cz;
  PetscReal    val;
  MatNullSpace nullspace;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Zero the incoming vector.
  PetscCall(VecZeroEntries(b));

  // Get the DM associated with the KSP.
  PetscCall(KSPGetDM(ksp, &dm));

  // Get an array equivalent to the RHS Vec.
  PetscCall(DMDAVecGetArray(dm, b, &rhs));

  // Get this processor's indices.
  PetscCall(DMDAGetCorners(dm, &i0, &j0, &k0, &ni, &nj, &nk));

  // Loop over grid points.
  for (k=k0; k<k0+nk; k++) {
    for (j=j0; j<j0+nj; j++) {
      for (i=i0; i<i0+ni; i++) {
        x = ((PetscReal)i + 0.5)*dx;
        y = ((PetscReal)j + 0.5)*dy;
        z = ((PetscReal)k + 0.5)*dz;
        Cx = PetscCosScalar(2*PETSC_PI * x);
        Cy = PetscCosScalar(2*PETSC_PI * y);
        Cz = PetscCosScalar(2*PETSC_PI * z);
        val = Cx * Cy * Cz;
        rhs[k][j][i] = n0 * val * scale;
      }
    }
  }

  // Restore the borrowed RHS array.
  PetscCall(DMDAVecRestoreArray(dm, b, &rhs));

  // Make the RHS vector consistent with the LHS operator.
  PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
  PetscCall(MatNullSpaceRemove(nullspace, b));
  PetscCall(MatNullSpaceDestroy(&nullspace));

  // Store the RHS vector in the problem context.
  ctx->potential.forcing = b;

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeFullRHS(KSP ksp, Vec b, void *user)
{
  Context          *ctx=(Context *)user;
  PetscReal         dx=ctx->grid.dx;
  PetscReal         dy=ctx->grid.dy;
  PetscReal         dz=ctx->grid.dz;
  PetscReal         scale=ctx->potential.scale;
  PetscReal         kappa=ctx->electrons.kappa;
  DM                fluidDM=ctx->fluidDM;
  Vec               moments;
  FluidNode      ***fluid;
  DM                dm;
  PetscReal      ***rhs;
  PetscInt          i0, j0, k0;
  PetscInt          ni, nj, nk;
  PetscInt          i, j, k;
  PetscReal         nijk, npjk, nmjk, nipk, nimk, nijp, nijm;
  PetscReal         hx, hy, hz;
  PetscReal         hxx, hyy, hzz;
  PetscReal         dndx, dndy;
  PetscReal         dGxdx, dGydy, dGzdz;
  PetscReal         d2ndxx, d2ndyy, d2ndzz;
  PetscReal         E0=ctx->plasma.E0;
  PetscReal         cth;
  PetscReal         cG;
  PetscReal         Eterm=0.0, Pterm=0.0, Gterm=0.0;
  MatNullSpace      nullspace;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Precompute scale factors for each term.
  cth = ctx->gammaT * KB * ctx->electrons.T / Q;
  cG =(1 + kappa*kappa) * ctx->electrons.m * ctx->electrons.nu / Q;

  // Precompute differential scale factors.
  hx = 1.0 / (2.0*dx);
  hy = 1.0 / (2.0*dy);
  hz = 1.0 / (2.0*dz);
  hxx = 1.0 / (dx*dx);
  hyy = 1.0 / (dy*dy);
  hzz = 1.0 / (dz*dz);

  // Get density and flux fluids.
  PetscCall(DMGetLocalVector(fluidDM, &moments));
  PetscCall(DMGlobalToLocalBegin(fluidDM, ctx->moments, INSERT_VALUES, moments));
  PetscCall(DMGlobalToLocalEnd(fluidDM, ctx->moments, INSERT_VALUES, moments));
  PetscCall(DMDAVecGetArray(fluidDM, moments, &fluid));

  // Zero the incoming vector.
  PetscCall(VecZeroEntries(b));

  // Get the DM associated with the KSP.
  PetscCall(KSPGetDM(ksp, &dm));

  // Get an fluid equivalent to the RHS Vec.
  PetscCall(DMDAVecGetArray(dm, b, &rhs));

  // Get this processor's indices.
  PetscCall(DMDAGetCorners(dm, &i0, &j0, &k0, &ni, &nj, &nk));

  // Loop over grid points.
  for (k=k0; k<k0+nk; k++) {
    for (j=j0; j<j0+nj; j++) {
      for (i=i0; i<i0+ni; i++) {

        /* Extract local density values for reuse. */
        nijk = fluid[k][j][i].n;
        npjk = fluid[k][j][i+1].n;
        nmjk = fluid[k][j][i-1].n;
        nipk = fluid[k][j+1][i].n;
        nimk = fluid[k][j-1][i].n;
        nijp = fluid[k+1][j][i].n;
        nijm = fluid[k-1][j][i].n;

        /* Compute local derivatives. */
        dndx = hx*(npjk - nmjk);
        dndy = hy*(nipk - nimk);
        d2ndxx = hxx*(npjk - 2*nijk + nmjk);
        d2ndyy = hyy*(nipk - 2*nijk + nimk);
        d2ndzz = hzz*(nijp - 2*nijk + nijm);
        dGxdx = hx*(fluid[k][j][i+1].Gx - fluid[k][j][i-1].Gx);
        dGydy = hy*(fluid[k][j+1][i].Gy - fluid[k][j-1][i].Gy);
        dGzdz = hz*(fluid[k+1][j][i].Gz - fluid[k-1][j][i].Gz);

        /* Assign the RHS value at (i, j, k). */
        Eterm = // div(n R E0)
          (dndy - kappa*dndx)*E0;
        Pterm = // (me vTe / e) div(R div(P))
          cth * (d2ndxx + d2ndyy + (1 + kappa*kappa)*d2ndzz);
        Gterm = // (1+kappa^2) (me nue / e) div(flux)
          cG * (dGxdx + dGydy + dGzdz);
        rhs[k][j][i] = scale * (Eterm + Pterm + Gterm);
      }
    }
  }

  // Restore density and flux fluids.
  PetscCall(DMDAVecRestoreArray(fluidDM, moments, &fluid));
  PetscCall(DMRestoreLocalVector(fluidDM, &moments));

  // Restore the borrowed RHS array.
  PetscCall(DMDAVecRestoreArray(dm, b, &rhs));

  // Make the RHS vector consistent with the LHS operator.
  PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullspace));
  PetscCall(MatNullSpaceRemove(nullspace, b));
  PetscCall(MatNullSpaceDestroy(&nullspace));

  // Store the RHS vector in the problem context.
  ctx->potential.forcing = b;

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


