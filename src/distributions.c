#include <petsc.h>
#include "hybrid.h"
#include "random.h"
#include "distributions.h"


PetscErrorCode UniformDistributionFromCoordinates(Context *ctx)
{

  PetscReal min[NDIM], max[NDIM];
  PetscInt  npc=1;                // number of particles per cell
  PetscInt  npd[NDIM];
  DM        swarmDM=ctx->swarmDM;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  min[0] = ctx->grid.x0 + 0.0*ctx->grid.dx;
  min[1] = ctx->grid.y0 + 0.0*ctx->grid.dy;
  min[2] = ctx->grid.z0 + 0.0*ctx->grid.dz;
  max[0] = ctx->grid.x1 + 0.0*ctx->grid.dx;
  max[1] = ctx->grid.y1 + 0.0*ctx->grid.dy;
  max[2] = ctx->grid.z1 + 0.0*ctx->grid.dz;
  npd[0] = npc * ctx->grid.Nx;
  npd[1] = npc * ctx->grid.Ny;
  npd[2] = npc * ctx->grid.Nz;

  // Use a built-in PETSc routine for setting up a uniform distribution.
  PetscCall(DMSwarmSetPointsUniformCoordinates(swarmDM, min, max, npd, INSERT_VALUES));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode UniformDistributionCellCentered(Context *ctx)
{
  PetscInt     npc=1;
  PetscInt     npt;
  DM           swarmDM=ctx->swarmDM;
  DM           cellDM;
  PetscScalar *pos;
  PetscInt     i0, j0, k0;
  PetscInt     ni, nj, nk;
  PetscInt     i, j, k;
  PetscInt     ip;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get information about the discrete grid.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));
  PetscCall(DMDAGetCorners(cellDM, &i0, &j0, &k0, &ni, &nj, &nk));

  // Compute the total number of particles.
  npt = ni*nj*nk*npc;

  // Reset the local swarm size to avoid a seg fault when accessing the
  // coordinates array. Passing a negative value for the buffer forces the swarm
  // to use its existing buffer size.
  PetscCall(DMSwarmSetLocalSizes(swarmDM, npt, -1));

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Loop over grid cells.
  ip = 0;
  for (k=k0; k<nk; k++) {
    for (j=j0; j<nj; j++) {
      for (i=i0; i<ni; i++) {
        pos[ip*NDIM + 0] = (PetscReal)i + 0.5;
        pos[ip*NDIM + 1] = (PetscReal)j + 0.5;
        pos[ip*NDIM + 2] = (PetscReal)k + 0.5;
        ip++;
      }
    }
  }

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode UniformDistribution(Context *ctx)
{
  DM           fluidDM=ctx->fluidDM;
  DM           swarmDM=ctx->swarmDM;
  PetscScalar *coords;
  PetscInt     np, np_cell, ip;
  PetscInt     i0, j0, k0;
  PetscInt     ni, nj, nk, nc;
  PetscInt     i, j, k, idx;
  PetscInt     dim;
  PetscReal    r[NDIM];
  PetscReal    d[NDIM]={ctx->grid.dx, ctx->grid.dy, ctx->grid.dz};

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get information about the discrete grid.
  PetscCall(DMDAGetCorners(fluidDM, &i0, &j0, &k0, &ni, &nj, &nk));

  // Get the local number of ions.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Compute the number of ions per cell. Note that np_cell*nc will not in
  // general be equal to the input value of -Np, if given.
  nc = ni*nj*nk;
  np_cell = (PetscInt)(np / nc);

  // Reset the local swarm size to avoid a seg fault when accessing the
  // coordinates array. Passing a negative value for the buffer forces the swarm
  // to use its existing buffer size.
  PetscCall(DMSwarmSetLocalSizes(swarmDM, np_cell*nc, -1));

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Loop over cells and place an equal number of particles at the center of
  // each cell in all but the last row of each index. This will result in a
  // symmetric step-function distribution. It is not truly uniform.
  for (ip=0; ip<=np_cell; ip++) {
    for (i=i0; i<i0+ni-1; i++) {
      for (j=j0; j<j0+nj-1; j++) {
        for (k=k0; k<k0+nk-1; k++) {
          idx = (ip*nc + k + j*nk + i*nk*nj)*NDIM;
          r[0] = (PetscReal)(i + 1);
          r[1] = (PetscReal)(j + 1);
          r[2] = (PetscReal)(k + 1);
          for (dim=0; dim<NDIM; dim++) {
            coords[idx + dim] = d[dim]*(r[dim] - 0.5);
          }
        }
      }
    }
  }

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode SobolDistribution(Context *ctx)
{
  DM             swarmDM=ctx->swarmDM;
  DM             cellDM;
  PetscInt       seed=-1, ndim=NDIM;
  PetscReal     *coords;
  PetscInt       Np, np, ip, ic;
  PetscReal      r[NDIM];
  PetscReal      lmin[NDIM], lmax[NDIM];
  PetscReal      d[NDIM]={ctx->grid.dx, ctx->grid.dy, ctx->grid.dz};
  PetscReal      L[NDIM]={ctx->grid.Lx, ctx->grid.Ly, ctx->grid.Lz};
  PetscInt       dim;
  DMDALocalInfo  local;
  PetscReal     *pos, x, y, z;
  PetscReal      xmin, xmax, ymin, ymax, zmin, zmax;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the total number of particles in the swarm.
  PetscCall(DMSwarmGetSize(swarmDM, &Np));

  // Allocate a 1-D array for the global positions.
  PetscCall(PetscMalloc1(NDIM*Np, &pos));

  if (ctx->mpi.rank == 0) {

    // Initialize the psuedo-random number generator.
    PetscCall(Sobseq(&seed, r-1));

    // Generate a Sobol' sequence of global positions on rank 0.
    for (ip=0; ip<Np; ip++) {
      PetscCall(Sobseq(&ndim, r-1));
      for (dim=0; dim<NDIM; dim++) {
        pos[ip*NDIM + dim] = r[dim]*L[dim];
      }
    }

  }

  // Broadcast the global positions array.
  PetscCallMPI(MPI_Bcast(pos, NDIM*Np, MPIU_REAL, 0, PETSC_COMM_WORLD));

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Get the local number of particles.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Get the ion-swarm cell DM.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));

  // Get the index information for this processor.
  PetscCall(DMDAGetLocalInfo(cellDM, &local));
  lmin[0] = d[0]*local.xs;
  lmin[1] = d[1]*local.ys;
  lmin[2] = d[2]*local.zs;
  lmax[0] = d[0]*(local.xs + local.xm);
  lmax[1] = d[1]*(local.ys + local.ym);
  lmax[2] = d[2]*(local.zs + local.zm);

  // Loop over global positions and assign local positions for this rank.
  /*
    NOTE: I think we can get away with a single loop here (as opposed to one
    loop to count the number of local particles, followed by allocating the
    local array, then a second loop to assign local positions) because DMSwarm
    has already allocated space for the estimated number of local particles.
  */
  xmin = lmin[0];
  xmax = lmax[0];
  ymin = lmin[1];
  ymax = lmax[1];
  zmin = lmin[2];
  zmax = lmax[2];
  ic = 0;
  for (ip=0; ip<Np; ip++) {
    x = pos[ip*NDIM + 0];
    y = pos[ip*NDIM + 1];
    z = pos[ip*NDIM + 2];
    if ((xmin <= x) && (x < xmax) && (ymin <= y) && (y < ymax) && (zmin <= z) && (z < zmax)) {
      coords[ic*NDIM + 0] = x;
      coords[ic*NDIM + 1] = y;
      coords[ic*NDIM + 2] = z;
      ic++;
    }
  }

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Free the positions-array memory.
  PetscCall(PetscFree(pos));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode SinusoidalDistribution(PetscReal x, PetscReal y, PetscReal z, PetscReal *v, Context *ctx)
{
  PetscReal fx;

  PetscFunctionBeginUser;

  fx = PetscSinReal(2*PETSC_PI * x);
  *v = 1.0 + 0.25*fx;

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode Rejection(DistributionFunction density, Context *ctx)
{
  PetscRandom  random;
  DM           fluidDM=ctx->fluidDM;
  DM           swarmDM=ctx->swarmDM;
  PetscInt     np, ip;
  PetscInt     i0, ni, i;
  PetscInt     j0, nj, j;
  PetscInt     k0, nk, k;
  PetscReal    localMax=0.0;
  PetscScalar *coords;
  PetscReal    Lx=ctx->grid.Lx;
  PetscReal    Ly=ctx->grid.Ly;
  PetscReal    Lz=ctx->grid.Lz;
  PetscReal    x, y, z, v, w;
  PetscReal    r[NDIM];
  PetscInt     it=0;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get a representation of the ion coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Create a random number generator.
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &random));
  PetscCall(PetscRandomSetInterval(random, 0.0, 1.0));
  PetscCall(PetscRandomSetSeed(random, ctx->mpi.rank));
  PetscCall(PetscRandomSeed(random));

  // Compute local maximum density.
  PetscCall(DMDAGetCorners(fluidDM, &i0, &j0, &k0, &ni, &nj, &nk));
  for (i=i0; i<i0+ni; i++) {
    for (j=j0; j<j0+nj; j++) {
      for (k=k0; k<k0+nk; k++) {
        PetscCall(density(i, j, k, &w, ctx));
        localMax = PetscMax(localMax, w);
      }
    }
  }
  ctx->log.ranks("[%d] Local maximum density: %g\n", ctx->mpi.rank, localMax);

  // Get the local number of ions.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Loop over all local ions.
  ip = 0;
  while (ip < np) {
    PetscCall(PetscRandomGetValuesReal(random, NDIM, r));
    x = r[0] * Lx;
    y = r[1] * Ly;
    z = r[2] * Lz;
    PetscCall(density(x, y, z, &w, ctx));
    PetscCall(PetscRandomGetValueReal(random, &v));
    if (w > v * localMax) {
      coords[ip*NDIM + 0] = x;
      coords[ip*NDIM + 1] = y;
      coords[ip*NDIM + 2] = z;
      ip++;
    }
    it++;
  }

  // Echo rejection efficiency.
  ctx->log.world("\n");
  ctx->log.ranks("[%d] Rejection efficiency: %f\n", ctx->mpi.rank, (PetscReal)ip/it);

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Destroy the random-number generator.
  PetscCall(PetscRandomDestroy(&random));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


