#include <petsc.h>
#include "ecsprt.h"
#include "random.h"
#include "velocities.h"


/* Generate a normal distribution of velocities with zero mean and unit variance. */
PetscErrorCode NormalVelocities(Context *ctx)
{
  DM         swarmDM=ctx->swarmDM;
  PetscInt   np, ip;
  PetscReal *vel;
  PetscReal  dvx, dvy, dvz;
  long       seed=getseed(*ctx);

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the number of local ions.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Loop over ions and assign parameter values.
  for (ip=0; ip<np; ip++) {
    PetscCall(Gasdev(&seed, &dvx));
    PetscCall(Gasdev(&seed, &dvy));
    PetscCall(Gasdev(&seed, &dvz));
    vel[ip*NDIM + 0] = ctx->ions.vTx*dvx + ctx->ions.v0x;
    vel[ip*NDIM + 1] = ctx->ions.vTy*dvy + ctx->ions.v0y;
    vel[ip*NDIM + 2] = ctx->ions.vTz*dvz + ctx->ions.v0z;
  }

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}

