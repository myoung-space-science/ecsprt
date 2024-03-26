#include "ecsprt.h"


/* Apply boundary conditions to particle positions and velocities. */
PetscErrorCode ApplyBC(PetscInt ndim, PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  Context   *ctx=(Context *)opts;
  PetscInt   ip, dim;
  PetscReal  r[ndim];
  BCType     bc[3]={ctx->ions.xBC, ctx->ions.yBC, ctx->ions.zBC};
  PetscReal  L[3]={ctx->grid.Lx, ctx->grid.Ly, ctx->grid.Lz};
  PetscReal  r0[3]={ctx->grid.x0, ctx->grid.y0, ctx->grid.z0};
  PetscReal  r1[3]={ctx->grid.x1, ctx->grid.y1, ctx->grid.z1};

  PetscFunctionBeginUser;

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    for (dim=0; dim<ndim; dim++) {
      // Extract the current coordinate.
      r[dim] = pos[ip*ndim + dim];
      // Check periodic boundary conditions.
      if (bc[dim] == BC_PERIODIC) {
        if (r[dim] < r0[dim]) {r[dim] += L[dim];}
        if (r[dim] >= r1[dim]) {r[dim] -= L[dim];}
      }
      // Copy the updated coordinate.
      pos[ip*ndim + dim] = r[dim];
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode Apply2DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  PetscFunctionBeginUser;
  PetscCall(ApplyBC(2, np, pos, vel, opts));
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode Apply3DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  PetscFunctionBeginUser;
  PetscCall(ApplyBC(3, np, pos, vel, opts));
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


