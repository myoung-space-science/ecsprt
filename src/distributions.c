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
  ECHO_FUNCTION_ENTER;

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

  ECHO_FUNCTION_EXIT;
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
  ECHO_FUNCTION_ENTER;

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

  ECHO_FUNCTION_EXIT;
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
  ECHO_FUNCTION_ENTER;

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

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode SobolDistribution(Context *ctx)
{
  DM          swarmDM=ctx->swarmDM;
  PetscInt    seed=-1, ndim=NDIM;
  PetscReal   *coords;
  PetscInt    np, ip;
  PetscReal   r[NDIM];
  PetscReal   L[NDIM]={ctx->grid.Lx, ctx->grid.Ly, ctx->grid.Lz};
  PetscInt    dim;

  PetscFunctionBeginUser;
  ECHO_FUNCTION_ENTER;

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Initialize the psuedo-random number generator.
  PetscCall(Sobseq(&seed, r-1));

  // Get the local number of particles.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  for (ip=0; ip<np; ip++) {
    PetscCall(Sobseq(&ndim, r-1));
    for (dim=0; dim<NDIM; dim++) {
      coords[ip*NDIM + dim] = r[dim]*L[dim];
    }
  }

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  ECHO_FUNCTION_EXIT;
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
  ECHO_FUNCTION_ENTER;

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
  PRINT_RANKS("[%d] Local maximum density: %g\n", ctx->mpi.rank, localMax);

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
  NEWLINE;
  PRINT_RANKS("[%d] Rejection efficiency: %f\n", ctx->mpi.rank, (PetscReal)ip/it);

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Destroy the random-number generator.
  PetscCall(PetscRandomDestroy(&random));

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


