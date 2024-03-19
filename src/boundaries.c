#include "ecsprt.h"


PetscErrorCode Apply2DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  Context   *ctx=(Context *)opts;
  PetscInt   ip, ndim=2;
  PetscReal  x, y;
  PetscReal  Lx=ctx->grid.Lx;
  PetscReal  Ly=ctx->grid.Ly;
  PetscReal  x0=ctx->grid.x0;
  PetscReal  y0=ctx->grid.y0;
  PetscReal  x1=ctx->grid.x1;
  PetscReal  y1=ctx->grid.y1;

  PetscFunctionBeginUser;

  // TODO: Implement other BC in ion loop. Some of those BC will require
  // modifying the velocity (e.g., vx = -vx for reflection at the x boundary, or
  // vx = vx0 for injection/advection along the x axis).

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    // Extract current positions.
    x = pos[ip*ndim + 0];
    y = pos[ip*ndim + 1];
    // Update the x position.
    if (ctx->ions.xBC == BC_PERIODIC) {
      if (x < x0) {x += Lx;}
      if (x >= x1) {x -= Lx;}
    }
    // Update the y position.
    if (ctx->ions.yBC == BC_PERIODIC) {
      if (y < y0) {y += Ly;}
      if (y >= y1) {y -= Ly;}
    }
    // Copy new positions.
    pos[ip*ndim + 0] = x;
    pos[ip*ndim + 1] = y;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode Apply3DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  Context   *ctx=(Context *)opts;
  PetscInt  ip, ndim=3;
  PetscReal x, y, z;
  PetscReal Lx=ctx->grid.Lx;
  PetscReal Ly=ctx->grid.Ly;
  PetscReal Lz=ctx->grid.Lz;
  PetscReal x0=ctx->grid.x0;
  PetscReal y0=ctx->grid.y0;
  PetscReal z0=ctx->grid.z0;
  PetscReal x1=ctx->grid.x1;
  PetscReal y1=ctx->grid.y1;
  PetscReal z1=ctx->grid.z1;

  PetscFunctionBeginUser;

  // TODO: Implement other BC in ion loop. Some of those BC will require
  // modifying the velocity (e.g., vx = -vx for reflection at the x boundary, or
  // vx = vx0 for injection/advection along the x axis).

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    // Extract current positions.
    x = pos[ip*ndim + 0];
    y = pos[ip*ndim + 1];
    z = pos[ip*ndim + 2];
    // Update the x position.
    if (ctx->ions.xBC == BC_PERIODIC) {
      if (x < x0) {x += Lx;}
      if (x >= x1) {x -= Lx;}
    }
    // Update the y position.
    if (ctx->ions.yBC == BC_PERIODIC) {
      if (y < y0) {y += Ly;}
      if (y >= y1) {y -= Ly;}
    }
    // Update the z position.
    if (ctx->ions.zBC == BC_PERIODIC) {
      if (z < z0) {z += Lz;}
      if (z >= z1) {z -= Lz;}
    }
    // Copy new positions.
    pos[ip*ndim + 0] = x;
    pos[ip*ndim + 1] = y;
    pos[ip*ndim + 2] = z;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Enforce boundary conditions and migrate particles.

This routine should be the only place where the code calls `DMSwarmMigrate`.
*/
PetscErrorCode ApplyBCAndMigrate(Context *ctx)
{
  DM         swarmDM=ctx->swarmDM;
  PetscReal *pos, *vel;
  PetscInt   np, Np0, Np1;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get an array representation of the ion positions.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Get the number of particles on this rank.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Apply boundary conditions.
  PetscCall(ctx->ions.applyBC(np, pos, vel, ctx));

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Echo pre-migration sizes.
  PetscCall(DMSwarmGetSize(swarmDM, &Np0));
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));
  ctx->log.world("\n");
  ctx->log.ranks("[%d] Local # of ions before migration: %d\n", ctx->mpi.rank, np);
  ctx->log.world("   Global # of ions before migration: %d\n", Np0);

  // Migrate ions among processes.
  PetscCall(DMSwarmMigrate(swarmDM, PETSC_TRUE));

  // Echo post-migration sizes.
  PetscCall(DMSwarmGetSize(swarmDM, &Np1));
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));
  ctx->log.world("\n");
  ctx->log.ranks("[%d] Local # of ions after migration: %d\n", ctx->mpi.rank, np);
  ctx->log.world("   Global # of ions after migration: %d\n", Np1);

  if ((ctx->ions.xBC == BC_PERIODIC) && (ctx->ions.yBC == BC_PERIODIC) && (ctx->ions.zBC == BC_PERIODIC)) {
    if (Np1 != Np0) {
      ctx->log.status("\n");
      ctx->log.status("ERROR: Total number of ions has changed (%+d).\n", Np1-Np0);
      ctx->log.status("       This should not happen with fully periodic boundary conditions.\n");
      ctx->log.status("\n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


