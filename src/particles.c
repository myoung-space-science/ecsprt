#include <petsc.h>
#include "ecsprt.h"
#include "positions.h"
#include "velocities.h"
#include "random.h"
#include "vector-math.h"
#include "particles.h"
#include "pic.h"
#include "boundaries.h"
#include "push.h"
#include "collisions.h"

/* Names of supported density functions. */
const char *PDistTypes[] ={
  "flat-sobol", "flat-reverse", "flat-normal", "uniform", "uniform-coordinates", "uniform-centered", "sinusoidal", "gaussian", "PDistType", "PDIST_", NULL
};

/* Names of supported initial velocity distributions. */
const char *VDistTypes[] = {
  "normal", "VDistType", "VDIST_", NULL
};

/* Names of supported boundary conditions. */
const char *BCTypes[] = {
  "periodic", "injection-reflection", "BC_", NULL
};


/* Compute the initial ion positions. */
PetscErrorCode InitializePositions(PDistType PDistType, Context *ctx)
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
    case PDIST_FLAT_NORMAL:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Not implemented: %s density", PDistTypes[PDIST_FLAT_NORMAL]);
      break;
    case PDIST_FLAT_REVERSE:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Not implemented: %s density", PDistTypes[PDIST_FLAT_REVERSE]);
      break;
    case PDIST_FLAT_SOBOL:
      PetscCall(SobolDistribution(ctx));
      break;
    case PDIST_UNIFORM:
      PetscCall(UniformDistribution(ctx));
      break;
    case PDIST_UNIFORM_COORDINATES:
      PetscCall(UniformDistributionFromCoordinates(ctx));
      break;
    case PDIST_UNIFORM_CENTERED:
      PetscCall(UniformDistributionCellCentered(ctx));
      break;
    case PDIST_SINUSOIDAL:
      PetscCall(Rejection(SinusoidalDistribution, ctx));
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


/* Compute the initial ion velocities. */
PetscErrorCode InitializeVelocities(VDistType VDistType, Context *ctx)
{
  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  switch (VDistType) {
    case VDIST_NORMAL:
      PetscCall(NormalVelocities(ctx));
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown initial velocity distribution");
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Update the ion velocities based on external forces and collisions. */
PetscErrorCode UpdateVelocities(PetscReal dt, Context *ctx)
{
  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Apply the Boris mover to integrate dv/dt = E + vxB.
  PetscCall(BorisMoverBz(dt, ctx));

  // Apply the appropriate collision algorithm.
  if (ctx->neutrals.m > 0.0) {
    PetscCall(ctx->ions.collide(dt, ctx));
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Update the ion positions according to $\frac{d\vec{r}}{dt} = \vec{v}$. */
PetscErrorCode UpdatePositions(PetscReal dt, Context *ctx)
{
  DM          swarmDM=ctx->swarmDM;
  PetscReal  *pos, *vel;
  PetscInt    ip, np;
  PetscReal   x, y, z;

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
    // Update the x position.
    x = pos[ip*NDIM + 0] + vel[ip*NDIM + 0]*dt;
    // Update the y position.
    y = pos[ip*NDIM + 1] + vel[ip*NDIM + 1]*dt;
    // Update the z position.
    z = pos[ip*NDIM + 2] + vel[ip*NDIM + 2]*dt;
    // Copy new positions.
    pos[ip*NDIM + 0] = x;
    pos[ip*NDIM + 1] = y;
    pos[ip*NDIM + 2] = z;
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


