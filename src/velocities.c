#include <petsc.h>
#include "ecsprt.h"
#include "random.h"
#include "velocities.h"
#include "vector-math.h"

/* Names of supported initial velocity distributions. */
const char *VDistTypes[] = {
  "normal", "VDistType", "VDIST_", NULL
};


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


/* Apply the 3-D vector Boris-mover algorithm to ion velocities.

See latter part of Birdsall & Langdon section 4-4.
*/
PetscErrorCode BorisMoverBB(PetscReal dt, Context *ctx)
{
  PetscReal   q=ctx->ions.q;
  PetscReal   m=ctx->ions.m;
  PetscReal   B[NDIM]={0.0, 0.0, ctx->plasma.B0};
  PetscInt    dim;
  PetscReal   dx=ctx->grid.dx, dy=ctx->grid.dy, dz=ctx->grid.dz;
  PetscReal   h[NDIM]={1.0/dx, 1.0/dy, 1.0/dz};
  PetscReal   t[NDIM], s[NDIM], t_dot_t;
  PetscReal   tscale, Escale[NDIM];
  DM          swarmDM=ctx->swarmDM;
  DM          phiDM=ctx->potential.dm;
  Vec         phiGlobal=ctx->potential.solution, phiLocal;
  PetscReal   ***phi;
  PetscInt    np, ip;
  PetscReal   *pos, *vel;
  PetscReal   x, y, z;
  PetscReal   E[NDIM]={0.0, 0.0, 0.0};
  PetscReal   vminus[NDIM], vprime[NDIM], vplus[NDIM];
  PetscReal   vminus_cross_t[NDIM], vprime_cross_s[NDIM];

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  /* Compute \vec{t} = \frac{q\vec{B}}{m}\frac{\Delta t}{2}. */
  tscale = 0.5 * (q/m) * dt;
  for (dim=0; dim<NDIM; dim++) {
    t[dim] = tscale * B[dim];
  }

  /* Compute \vec{s} = \frac{2\vec{t}}{1 + \vec{t}\cdot\vec{t}}. */
  PetscCall(DotProduct(t, t, &t_dot_t));
  for (dim=0; dim<NDIM; dim++) {
    s[dim] = 2.0 * t[dim] / (1 + t_dot_t);
  }

  /* Compute the electric-field scale factors.

  These account for the species constants as well as the 2nd-order
  finite-difference gradient scale factors.
  */
  for (dim=0; dim<NDIM; dim++) {
    Escale[dim] = -0.5 * h[dim] * tscale;
  }

  /* Get a local copy of phi with ghost cells. */
  PetscCall(DMGetLocalVector(phiDM, &phiLocal));
  PetscCall(DMGlobalToLocal(phiDM, phiGlobal, INSERT_VALUES, phiLocal));

  /* Get a temporary array representing the local electrostatic potential. */
  PetscCall(DMDAVecGetArray(phiDM, phiLocal, &phi));

  /* Get the number of local ions. */
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  /* Get an array representation of the ion positions. */
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  /* Get an array representation of the ion velocities. */
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  /* Loop over particles and interpolate E to grid points. */
  for (ip=0; ip<np; ip++) {
    /* Get the current particle's coordinates. */
    /* Normalize each coordinate to a fractional number of grid cells. */
    x = pos[ip*NDIM + 0] / dx;
    y = pos[ip*NDIM + 1] / dy;
    z = pos[ip*NDIM + 2] / dz;
    /* Compute the electric field due to this particle: \vec{E} = -\nabla\phi. */
    PetscCall(DifferenceVector(phi, x, y, z, ctx->grid, E));
    /* Compute \vec{v}^-. */
    for (dim=0; dim<NDIM; dim++) {
      vminus[0] = vel[ip*NDIM + dim] + Escale[dim]*E[dim];
    }
    /* Compute \vec{v}^- \times \vec{t}. */
    PetscCall(CrossProduct(vminus, t, vminus_cross_t));
    /* Compute \vec{v}^\prime = \vec{v}^- + \vec{v}^- \times \vec{t}. */
    for (dim=0; dim<NDIM; dim++) {
      vprime[dim] = vminus[dim] + vminus_cross_t[dim];
    }
    /* Compute \vec{v}^\prime \times \vec{s}. */
    PetscCall(CrossProduct(vprime, s, vprime_cross_s));
    /* Compute \vec{v}^+ = \vec{v}^- + \vec{v}^\prime \times \vec{s}. */
    for (dim=0; dim<NDIM; dim++) {
      vplus[dim] = vminus[dim] + vprime_cross_s[dim];
    }
    /* Assign new particle velocities. */
    for (dim=0; dim<NDIM; dim++){
      vel[ip*NDIM + dim] = vplus[dim] + Escale[dim]*E[dim];
    }
  }

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Restore the borrowed potential array and local vector.
  PetscCall(DMDAVecRestoreArray(phiDM, phiLocal, &phi));
  PetscCall(DMRestoreLocalVector(phiDM, &phiLocal));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Apply the 3-D scalar Boris-mover algorithm to ion velocities.

See beginning of Birdsall & Langdon section 4-4. This algorithm assumes that
the magnetic field is aligned with the z dimension.
*/
PetscErrorCode BorisMoverBz(PetscReal dt, Context *ctx)
{
  PetscReal    q=ctx->ions.q;
  PetscReal    m=ctx->ions.m;
  PetscReal    B=ctx->plasma.B0;
  PetscInt     dim;
  PetscReal    dx=ctx->grid.dx, dy=ctx->grid.dy, dz=ctx->grid.dz;
  PetscReal    h[NDIM]={1.0/dx, 1.0/dy, 1.0/dz};
  PetscReal    t, s;
  PetscReal    tscale, Escale[NDIM];
  DM           swarmDM=ctx->swarmDM;
  DM           phiDM=ctx->potential.dm;
  Vec          phiGlobal=ctx->potential.solution, phiLocal;
  PetscReal ***phi;
  PetscInt     np, ip;
  PetscReal   *pos, *vel;
  PetscReal    x, y, z;
  PetscReal    E[NDIM]={0.0, 0.0, 0.0};
  PetscReal    vxm, vym, vz, vxt, vxp, vyp;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  /* Compute t = \frac{qB}{m}\frac{\Delta t}{2}. */
  tscale = 0.5 * (q/m) * dt;
  t = tscale * B;

  /* Compute s = \frac{2t}{1 + t^2}. */
  s = 2.0 * t / (1 + t*t);

  /* Compute the electric-field scale factors.

  These account for the species constants as well as the 2nd-order
  finite-difference gradient scale factors.
  */
  for (dim=0; dim<NDIM; dim++) {
    Escale[dim] = -0.5 * h[dim] * tscale;
  }

  /* Get a local copy of phi with ghost cells. */
  PetscCall(DMGetLocalVector(phiDM, &phiLocal));
  PetscCall(DMGlobalToLocal(phiDM, phiGlobal, INSERT_VALUES, phiLocal));

  /* Get a temporary array representing the local electrostatic potential. */
  PetscCall(DMDAVecGetArray(phiDM, phiLocal, &phi));

  /* Get the number of local ions. */
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  /* Get an array representation of the ion positions. */
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  /* Get an array representation of the ion velocities. */
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  /* Loop over particles and interpolate E to grid points. */
  for (ip=0; ip<np; ip++) {
    /* Get the current particle's coordinates. */
    /* Normalize each coordinate to a fractional number of grid cells. */
    x = pos[ip*NDIM + 0] / dx;
    y = pos[ip*NDIM + 1] / dy;
    z = pos[ip*NDIM + 2] / dz;
    /* Compute the electric field due to this particle: \vec{E} = -\nabla\phi. */
    PetscCall(DifferenceVector(phi, x, y, z, ctx->grid, E));
    vxm = vel[ip*NDIM + 0] + Escale[0]*E[0];
    vym = vel[ip*NDIM + 1] + Escale[1]*E[1];
    vz  = vel[ip*NDIM + 2] + Escale[2]*E[2];
    vxt = vxm + vym*t;
    vyp = vym - vxt*s;
    vxp = vxt + vyp*t;
    /* Assign new particle velocities. */
    vel[ip*NDIM + 0] = vxp;
    vel[ip*NDIM + 1] = vyp;
    vel[ip*NDIM + 2] = vz;
  }

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Restore the borrowed potential array and local vector.
  PetscCall(DMDAVecRestoreArray(phiDM, phiLocal, &phi));
  PetscCall(DMRestoreLocalVector(phiDM, &phiLocal));

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


