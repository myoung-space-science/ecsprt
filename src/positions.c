#include <petsc.h>
#include "ecsprt.h"
#include "random.h"
#include "positions.h"
#include "boundaries.h"

/* Names of supported density functions. */
const char *PDistTypes[] ={
  "sobol", "reverse", "normal", "uniform", "sinusoidal", "gaussian", "PDistType", "PDIST_", NULL
};


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


PetscErrorCode NormalDistribution(PetscInt ndim, Context *ctx)
{
  DM             swarmDM=ctx->swarmDM;
  DM             cellDM;
  long           seed=getseed(*ctx);
  PetscReal     *coords;
  PetscInt       Np, np, ip, ic;
  PetscReal      s;
  PetscReal      L[3]={ctx->grid.Lx, ctx->grid.Ly, ctx->grid.Lz};
  PetscReal      dx=ctx->grid.dx;
  PetscReal      dy=ctx->grid.dy;
  PetscReal      dz=ctx->grid.dz;
  PetscInt       dim;
  DMDALocalInfo  local;
  PetscReal     *pos, x, y, z;
  PetscReal      xmin, xmax, ymin, ymax, zmin, zmax;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the total number of particles in the swarm.
  PetscCall(DMSwarmGetSize(swarmDM, &Np));

  // Allocate a 1-D array for the global positions.
  PetscCall(PetscMalloc1(ndim*Np, &pos));

  if (ctx->mpi.rank == 0) {

    // Generate a normal distribution of global positions on rank 0.
    for (ip=0; ip<Np; ip++) {
      for (dim=0; dim<ndim; dim++) {
        PetscCall(Ran3(&seed, &s));
        pos[ip*ndim + dim] = s*L[dim];
      }
    }

  }

  // Broadcast the global positions array.
  PetscCallMPI(MPI_Bcast(pos, ndim*Np, MPIU_REAL, 0, PETSC_COMM_WORLD));

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Get the local number of particles.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Get the ion-swarm cell DM.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));

  // Get the index information for this processor.
  PetscCall(DMDAGetLocalInfo(cellDM, &local));
  xmin = dx*local.xs;
  ymin = dy*local.ys;
  zmin = dz*local.zs;
  xmax = dx*(local.xs + local.xm);
  ymax = dy*(local.ys + local.ym);
  zmax = dz*(local.zs + local.zm);

  // Loop over global positions and assign local positions for this rank.
  /* NOTE: DMSwarm has already allocated space for the coordinates array. */
  ic = 0;
  if (ndim == 2) {
    for (ip=0; ip<Np; ip++) {
      x = pos[ip*ndim + 0];
      y = pos[ip*ndim + 1];
      if ((xmin <= x) && (x < xmax) && (ymin <= y) && (y < ymax)) {
        coords[ic*ndim + 0] = x;
        coords[ic*ndim + 1] = y;
        ic++;
      }
    }
  } else if (ndim == 3) {
    for (ip=0; ip<Np; ip++) {
      x = pos[ip*ndim + 0];
      y = pos[ip*ndim + 1];
      z = pos[ip*ndim + 2];
      if ((xmin <= x) && (x < xmax) && (ymin <= y) && (y < ymax) && (zmin <= z) && (z < zmax)) {
        coords[ic*ndim + 0] = x;
        coords[ic*ndim + 1] = y;
        coords[ic*ndim + 2] = z;
        ic++;
      }
    }
  }

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Free the positions-array memory.
  PetscCall(PetscFree(pos));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode SobolDistribution(PetscInt ndim, Context *ctx)
{
  DM             swarmDM=ctx->swarmDM;
  DM             cellDM;
  PetscInt       seed=-1;
  PetscReal     *coords;
  PetscInt       Np, ip, ic;
  PetscReal      s[ndim];
  PetscReal      L[3]={ctx->grid.Lx, ctx->grid.Ly, ctx->grid.Lz};
  PetscReal      dx=ctx->grid.dx;
  PetscReal      dy=ctx->grid.dy;
  PetscReal      dz=ctx->grid.dz;
  PetscInt       dim;
  PetscInt       ntmp;
  DMDALocalInfo  local;
  PetscReal     *pos;
  PetscReal      r[3];
  PetscReal      lmin[3], lmax[3];

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the total number of particles in the swarm.
  PetscCall(DMSwarmGetSize(swarmDM, &Np));

  // Allocate a 1-D array for the global positions.
  PetscCall(PetscMalloc1(ndim*Np, &pos));

  // Generate a Sobol' sequence of global positions on rank 0.
  if (ctx->mpi.rank == 0) {

    PetscCall(Sobseq(&seed, s-1));
    for (ip=0; ip<Np; ip++) {
      PetscCall(Sobseq(&ndim, s-1));
      for (dim=0; dim<ndim; dim++) {
        pos[ip*ndim + dim] = s[dim]*L[dim];
      }
    }

  }

  // Broadcast the global positions array.
  PetscCallMPI(MPI_Bcast(pos, ndim*Np, MPIU_REAL, 0, PETSC_COMM_WORLD));

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Get the ion-swarm cell DM.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));

  // Get the index information for this processor.
  PetscCall(DMDAGetLocalInfo(cellDM, &local));
  lmin[0] = dx*local.xs;
  lmin[1] = dy*local.ys;
  lmin[2] = dz*local.zs;
  lmax[0] = dx*(local.xs + local.xm);
  lmax[1] = dy*(local.ys + local.ym);
  lmax[2] = dz*(local.zs + local.zm);

  // Loop over global positions and assign local positions for this rank.
  /* NOTE: DMSwarm has already allocated space for the coordinates array. */
  ic = 0;
  for (ip=0; ip<Np; ip++) {
    /* Store the current particle's global position. */
    for (dim=0; dim<ndim; dim++) {
      r[dim] = pos[ip*ndim + dim];
    }
    /* Test whether this particle's global position is within the local box. */
    ntmp = 0;
    for (dim=0; dim<ndim; dim++) {
      if ((lmin[dim] <= r[dim]) && (r[dim] < lmax[dim])) {
        ntmp++;
      }
    }
    /* If so, assign this particle's position to the set of local coordinates. */
    if (ntmp == ndim) {
      for (dim=0; dim<ndim; dim++) {
        coords[ic*ndim + dim] = r[dim];
      }
      ic++;
    }
  }

  // Restore the coordinates array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Reset the number of local particles.
  PetscCall(DMSwarmSetLocalSizes(swarmDM, ic, -1));

  // Free the positions-array memory.
  PetscCall(PetscFree(pos));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode SinusoidalDistribution(PetscReal x, PetscReal y, PetscReal z, PetscReal *v, Context *ctx)
{
  PetscReal fx;

  PetscFunctionBeginUser;

  // TODO: Define runtime parameters
  // - x, y, z amplitudes
  // - x, y, z periods
  // - x, y, z shifts
  fx = PetscSinReal(2*PETSC_PI * x);
  *v = 1.0 + 0.25*fx;

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode Rejection(DistributionFunction density, Context *ctx)
{
  DM             swarmDM=ctx->swarmDM;
  DM             cellDM;
  PetscReal     *coords;
  PetscInt       Np, ip, ic;
  PetscRandom    random;
  PetscInt       i, j, k;
  PetscInt       Nx=ctx->grid.Nx;
  PetscInt       Ny=ctx->grid.Ny;
  PetscInt       Nz=ctx->grid.Nz;
  PetscInt       Lx=ctx->grid.Lx;
  PetscInt       Ly=ctx->grid.Ly;
  PetscInt       Lz=ctx->grid.Lz;
  PetscReal      dx=ctx->grid.dx;
  PetscReal      dy=ctx->grid.dy;
  PetscReal      dz=ctx->grid.dz;
  PetscReal      maxVal=0.0, normVal;
  PetscReal      r[NDIM];
  PetscInt       it;
  DMDALocalInfo  local;
  PetscReal     *pos, x, y, z, w, v;
  PetscReal      xmin, xmax, ymin, ymax, zmin, zmax;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the total number of particles in the swarm.
  PetscCall(DMSwarmGetSize(swarmDM, &Np));

  // Allocate a 1-D array for the global positions.
  PetscCall(PetscMalloc1(NDIM*Np, &pos));

  if (ctx->mpi.rank == 0) {

    // Create a random number generator.
    PetscCall(PetscRandomCreate(PETSC_COMM_SELF, &random));
    PetscCall(PetscRandomSetInterval(random, 0.0, 1.0));
    PetscCall(PetscRandomSetSeed(random, 7));
    PetscCall(PetscRandomSeed(random));

    // Compute the maximum global density.
    for (i=0; i<Nx; i++) {
      for (j=0; j<Ny; j++) {
        for (k=0; k<Nz; k++) {
          x = (PetscReal)i / Nx;
          y = (PetscReal)j / Ny;
          z = (PetscReal)k / Nz;
          PetscCall(density(x, y, z, &w, ctx));
          maxVal = PetscMax(maxVal, w);
        }
      }
    }

    // Compute the global positions via rejection.
    normVal = 1.0 / maxVal;
    it = 0;
    while (ip < Np) {
      PetscCall(PetscRandomGetValuesReal(random, NDIM, r));
      x = r[0];
      y = r[1];
      z = r[2];
      PetscCall(density(x, y, z, &w, ctx));
      PetscCall(PetscRandomGetValueReal(random, &v));
      if (w*normVal > v) {
        pos[ip*NDIM + 0] = x*Lx;
        pos[ip*NDIM + 1] = y*Ly;
        pos[ip*NDIM + 2] = z*Lz;
        ip++;
      }
      it++;
    }

    // Destroy the random-number generator.
    PetscCall(PetscRandomDestroy(&random));

  }

  // Broadcast the global positions array.
  PetscCallMPI(MPI_Bcast(pos, NDIM*Np, MPIU_REAL, 0, PETSC_COMM_WORLD));

  // Get a representation of the particle coordinates.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&coords));

  // Get the ion-swarm cell DM.
  PetscCall(DMSwarmGetCellDM(swarmDM, &cellDM));

  // Get the index information for this processor.
  PetscCall(DMDAGetLocalInfo(cellDM, &local));
  xmin = dx*local.xs;
  xmax = dx*(local.xs + local.xm);
  ymin = dy*local.ys;
  ymax = dy*(local.ys + local.ym);
  zmin = dz*local.zs;
  zmax = dz*(local.zs + local.zm);

  // Loop over global positions and assign local positions for this rank.
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

  // Reset the number of local particles.
  PetscCall(DMSwarmSetLocalSizes(swarmDM, ic, -1));

  // Free the positions-array memory.
  PetscCall(PetscFree(pos));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the initial ion positions. */
PetscErrorCode InitializePositions(PetscInt ndim, PDistType PDistType, Context *ctx)
{
  DM         swarmDM=ctx->swarmDM;
  PetscInt   np, Np;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Echo sizes.
  PetscCall(DMSwarmGetSize(swarmDM, &Np));
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));
  ctx->log.world("\n");
  ctx->log.ranks("[%d] Local # of ions before placement: %d\n", ctx->mpi.rank, np);
  ctx->log.world("   Global # of ions before placement: %d\n", Np);

  // Initialize coordinates in the ions DM.
  switch(PDistType) {
    case PDIST_NORMAL:
      PetscCall(NormalDistribution(ndim, ctx));
      break;
    case PDIST_REVERSE:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Not implemented: %s density", PDistTypes[PDIST_REVERSE]);
      break;
    case PDIST_SOBOL:
      PetscCall(SobolDistribution(ndim, ctx));
      break;
    case PDIST_UNIFORM:
      PetscCall(UniformCoordinates(ndim, ctx));
      break;
    case PDIST_SINUSOIDAL:
      PetscCall(Rejection(ndim, SinusoidalDistribution, ctx));
      break;
    case PDIST_GAUSSIAN:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Not implemented: %s density", PDistTypes[PDIST_GAUSSIAN]);
      break;
  }

  // Apply BC and migrate particles.
  PetscCall(ApplyBCAndMigrate(ctx));

  // Update the parameter context.
  PetscCall(DMSwarmGetSize(swarmDM, &ctx->plasma.Np));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Update the ion positions according to $\frac{d\vec{r}}{dt} = \vec{v}$. */
PetscErrorCode UpdatePositions(PetscInt ndim, PetscReal dt, Context *ctx)
{
  DM          swarmDM=ctx->swarmDM;
  PetscReal  *pos, *vel;
  PetscInt    ip, np, dim;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get an array representation of the ion positions.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Get the number of particles on this rank.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    for (dim=0; dim<ndim; dim++) {
      pos[ip*ndim + dim] += vel[ip*ndim + dim]*dt;
    }
  }

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Update the swarm.
  PetscCall(ApplyBCAndMigrate(ctx));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


