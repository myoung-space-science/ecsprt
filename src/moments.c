#include "ecsprt.h"


/* Compute moments of a 2-D ion distribution. */
PetscErrorCode Collect2DFluidMoments(void *opts)
{
  Context     *ctx=(Context *)opts;
  PetscInt     ndim=2;
  DM           swarmDM=ctx->swarmDM;
  DM           cellDM;
  Vec          moments, global;
  PetscReal ***array;
  PetscReal    x, y, dx, dy;
  PetscReal   *pos;
  PetscReal   *vel;
  PetscInt     ip, np;
  PetscInt     ixl, ixh, iyl, iyh;
  PetscReal    wxl, wxh, wyl, wyh;
  PetscReal    hh, lh, hl, ll;
  PetscInt     dim, dof;
  PetscReal    w;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the ion-swarm cell DM.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));

  // Get density and flux arrays.
  PetscCall(DMGetLocalVector(cellDM, &moments));
  PetscCall(VecZeroEntries(moments));
  PetscCall(DMDAVecGetArrayDOF(cellDM, moments, &array));

  // Get an array representation of the ion positions.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Get the number of ions on this rank.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Extract cell widths for reuse.
  dx = ctx->grid.dx;
  dy = ctx->grid.dy;

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    // Normalize each coordinate to a fractional number of grid cells.
    x = pos[ip*ndim + 0] / dx;
    y = pos[ip*ndim + 1] / dy;

    // TODO: This needs to incorporate BC (at least for certain BC), in addition
    // to whatever action `ApplyBCAndMigrate` takes.

    // Compute the x-dimension neighbors and corresponding weights.
    ixl = (PetscInt)x;
    ixh = ixl+1;
    wxh = x - (PetscReal)ixl;
    wxl = 1.0 - wxh;
    // Compute the y-dimension neighbors and corresponding weights.
    iyl = (PetscInt)y;
    iyh = iyl+1;
    wyh = y - (PetscReal)iyl;
    wyl = 1.0 - wyh;
    // Compute the weight of each nearby grid point.
    hh = wyh*wxh;
    lh = wyl*wxh;
    hl = wyh*wxl;
    ll = wyl*wxl;
    // Assign density values (zeroth moment).
    array[iyh][ixh][0] += hh;
    array[iyl][ixh][0] += lh;
    array[iyh][ixl][0] += hl;
    array[iyl][ixl][0] += ll;
    // Assign flux values (first moments wrt velocity).
    for (dim=0, dof=1; dim<ndim; dim++, dof++) {
      w = vel[ip*ndim + dim];
      array[iyh][ixh][dof] += w*hh;
      array[iyl][ixh][dof] += w*lh;
      array[iyh][ixl][dof] += w*hl;
      array[iyl][ixl][dof] += w*ll;
    }
  }

  PetscCall(DMGetGlobalVector(cellDM, &global));
  PetscCall(VecZeroEntries(global));
  PetscCall(DMLocalToGlobal(cellDM, moments, ADD_VALUES, global));
  PetscCall(VecCopy(global, ctx->moments));
  PetscCall(DMRestoreGlobalVector(cellDM, &global));

  // Restore density and flux arrays.
  PetscCall(DMDAVecRestoreArrayDOF(cellDM, moments, &array));
  PetscCall(DMRestoreLocalVector(cellDM, &moments));

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute moments of a 3-D ion distribution. */
PetscErrorCode Collect3DFluidMoments(void *opts)
{
  Context      *ctx=(Context *)opts;
  PetscInt      ndim=3;
  DM            swarmDM=ctx->swarmDM;
  DM            cellDM;
  Vec           moments, global;
  PetscReal ****array;
  PetscReal     x, y, z, dx, dy, dz;
  PetscReal    *pos;
  PetscReal    *vel;
  PetscInt      ip, np;
  PetscInt      ixl, ixh, iyl, iyh, izl, izh;
  PetscReal     wxl, wxh, wyl, wyh, wzl, wzh;
  PetscReal     hhh, lhh, hlh, llh, hhl, lhl, hll, lll;
  PetscInt      dim, dof;
  PetscReal     w;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the ion-swarm cell DM.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));

  // Get density and flux arrays.
  PetscCall(DMGetLocalVector(cellDM, &moments));
  PetscCall(VecZeroEntries(moments));
  PetscCall(DMDAVecGetArrayDOF(cellDM, moments, &array));

  // Get an array representation of the ion positions.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Get the number of ions on this rank.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Extract cell widths for reuse.
  dx = ctx->grid.dx;
  dy = ctx->grid.dy;
  dz = ctx->grid.dz;

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    // Normalize each coordinate to a fractional number of grid cells.
    x = pos[ip*ndim + 0] / dx;
    y = pos[ip*ndim + 1] / dy;
    z = pos[ip*ndim + 2] / dz;

    // TODO: This needs to incorporate BC (at least for certain BC), in addition
    // to whatever action `ApplyBCAndMigrate` takes.

    // Compute the x-dimension neighbors and corresponding weights.
    ixl = (PetscInt)x;
    ixh = ixl+1;
    wxh = x - (PetscReal)ixl;
    wxl = 1.0 - wxh;
    // Compute the y-dimension neighbors and corresponding weights.
    iyl = (PetscInt)y;
    iyh = iyl+1;
    wyh = y - (PetscReal)iyl;
    wyl = 1.0 - wyh;
    // Compute the z-dimension neighbors and corresponding weights.
    izl = (PetscInt)z;
    izh = izl+1;
    wzh = z - (PetscReal)izl;
    wzl = 1.0 - wzh;
    // Compute the weight of each nearby grid point.
    hhh = wzh*wyh*wxh;
    lhh = wzl*wyh*wxh;
    hlh = wzh*wyl*wxh;
    llh = wzl*wyl*wxh;
    hhl = wzh*wyh*wxl;
    lhl = wzl*wyh*wxl;
    hll = wzh*wyl*wxl;
    lll = wzl*wyl*wxl;
    // Assign density values (zeroth moment).
    array[izh][iyh][ixh][0] += hhh;
    array[izl][iyh][ixh][0] += lhh;
    array[izh][iyl][ixh][0] += hlh;
    array[izl][iyl][ixh][0] += llh;
    array[izh][iyh][ixl][0] += hhl;
    array[izl][iyh][ixl][0] += lhl;
    array[izh][iyl][ixl][0] += hll;
    array[izl][iyl][ixl][0] += lll;
    // Assign flux values (first moments wrt velocity).
    for (dim=0, dof=1; dim<ndim; dim++, dof++) {
      w = vel[ip*ndim + dim];
      array[izh][iyh][ixh][dof] += w*hhh;
      array[izl][iyh][ixh][dof] += w*lhh;
      array[izh][iyl][ixh][dof] += w*hlh;
      array[izl][iyl][ixh][dof] += w*llh;
      array[izh][iyh][ixl][dof] += w*hhl;
      array[izl][iyh][ixl][dof] += w*lhl;
      array[izh][iyl][ixl][dof] += w*hll;
      array[izl][iyl][ixl][dof] += w*lll;
    }
  }

  PetscCall(DMGetGlobalVector(cellDM, &global));
  PetscCall(VecZeroEntries(global));
  PetscCall(DMLocalToGlobal(cellDM, moments, ADD_VALUES, global));
  PetscCall(VecCopy(global, ctx->moments));
  PetscCall(DMRestoreGlobalVector(cellDM, &global));

  // Restore density and flux arrays.
  PetscCall(DMDAVecRestoreArrayDOF(cellDM, moments, &array));
  PetscCall(DMRestoreLocalVector(cellDM, &moments));

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


