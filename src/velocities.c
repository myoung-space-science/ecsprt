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
PetscErrorCode NormalVelocities(PetscInt ndim, Context *ctx)
{
  DM         swarmDM=ctx->swarmDM;
  PetscInt   np, ip;
  PetscReal *vel;
  PetscReal  vT[3]={ctx->ions.vTx, ctx->ions.vTy, ctx->ions.vTz};
  PetscReal  v0[3]={ctx->ions.v0x, ctx->ions.v0y, ctx->ions.v0z};
  PetscReal  dv;
  long       seed=getseed(*ctx);
  PetscInt   dim;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get the number of local ions.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Loop over ions and assign parameter values.
  for (ip=0; ip<np; ip++) {
    for (dim=0; dim<ndim; dim++) {
      PetscCall(Gasdev(&seed, &dv));
      vel[ip*ndim + dim] = vT[dim]*dv + v0[dim];
    }
  }

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the initial ion velocities. */
PetscErrorCode InitializeVelocities(PetscInt ndim, VDistType VDistType, Context *ctx)
{
  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  switch (VDistType) {
    case VDIST_NORMAL:
      PetscCall(NormalVelocities(ndim, ctx));
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown initial velocity distribution");
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Apply the 3-D vector Boris-mover algorithm to ion velocities.

See latter part of Birdsall & Langdon section 4-4. Note that there is no 2-D
version of this routine.
*/
PetscErrorCode BorisMoverBB(PetscReal dt, Context *ctx)
{
  PetscReal   q=ctx->ions.q;
  PetscReal   m=ctx->ions.m;
  PetscReal   B[3]={0.0, 0.0, ctx->plasma.B0};
  PetscInt    dim;
  PetscReal   dx=ctx->grid.dx, dy=ctx->grid.dy, dz=ctx->grid.dz;
  PetscReal   h[3]={1.0/dx, 1.0/dy, 1.0/dz};
  PetscReal   t[3], s[3], t_dot_t;
  PetscReal   tscale, Escale[3];
  DM          swarmDM=ctx->swarmDM;
  DM          phiDM=ctx->potential.dm;
  Vec         phiGlobal=ctx->potential.solution, phiLocal;
  PetscReal   ***phi;
  PetscInt    np, ip;
  PetscReal   *pos, *vel;
  PetscReal   x, y, z;
  PetscReal   E[3]={0.0, 0.0, 0.0};
  PetscReal   vminus[3], vprime[3], vplus[3];
  PetscReal   vminus_cross_t[3], vprime_cross_s[3];

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  /* Compute \vec{t} = \frac{q\vec{B}}{m}\frac{\Delta t}{2}. */
  tscale = 0.5 * (q/m) * dt;
  for (dim=0; dim<3; dim++) {
    t[dim] = tscale * B[dim];
  }

  /* Compute \vec{s} = \frac{2\vec{t}}{1 + \vec{t}\cdot\vec{t}}. */
  PetscCall(DotProduct(t, t, &t_dot_t));
  for (dim=0; dim<3; dim++) {
    s[dim] = 2.0 * t[dim] / (1 + t_dot_t);
  }

  /* Compute the electric-field scale factors.

  These account for the species constants as well as the 2nd-order
  finite-difference gradient scale factors.
  */
  for (dim=0; dim<3; dim++) {
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
    x = pos[ip*3 + 0] / dx;
    y = pos[ip*3 + 1] / dy;
    z = pos[ip*3 + 2] / dz;
    /* Compute the electric field due to this particle: \vec{E} = -\nabla\phi. */
    PetscCall(DifferenceVector3D(phi, x, y, z, ctx->grid, E));
    /* Compute \vec{v}^-. */
    for (dim=0; dim<3; dim++) {
      vminus[0] = vel[ip*3 + dim] + Escale[dim]*E[dim];
    }
    /* Compute \vec{v}^- \times \vec{t}. */
    PetscCall(CrossProduct(vminus, t, vminus_cross_t));
    /* Compute \vec{v}^\prime = \vec{v}^- + \vec{v}^- \times \vec{t}. */
    for (dim=0; dim<3; dim++) {
      vprime[dim] = vminus[dim] + vminus_cross_t[dim];
    }
    /* Compute \vec{v}^\prime \times \vec{s}. */
    PetscCall(CrossProduct(vprime, s, vprime_cross_s));
    /* Compute \vec{v}^+ = \vec{v}^- + \vec{v}^\prime \times \vec{s}. */
    for (dim=0; dim<3; dim++) {
      vplus[dim] = vminus[dim] + vprime_cross_s[dim];
    }
    /* Assign new particle velocities. */
    for (dim=0; dim<3; dim++){
      vel[ip*3 + dim] = vplus[dim] + Escale[dim]*E[dim];
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


/* Apply the 2-D scalar Boris-mover algorithm to ion velocities.

See beginning of Birdsall & Langdon section 4-4. This algorithm assumes that
the magnetic field is aligned with the z dimension.
*/
PetscErrorCode BorisMoverBz2D(PetscReal dt, void *opts)
{
  Context    *ctx=(Context *)opts;
  PetscReal   q=ctx->ions.q;
  PetscReal   m=ctx->ions.m;
  PetscReal   B=ctx->plasma.B0;
  PetscInt    dim;
  PetscReal   dx=ctx->grid.dx, dy=ctx->grid.dy;
  PetscReal   h[2]={1.0/dx, 1.0/dy};
  PetscReal   t, s;
  PetscReal   tscale, Escale[2];
  DM          swarmDM=ctx->swarmDM;
  DM          phiDM=ctx->potential.dm;
  Vec         phiGlobal=ctx->potential.solution, phiLocal;
  PetscReal **phi;
  PetscInt    np, ip;
  PetscReal  *pos, *vel;
  PetscReal   x, y;
  PetscReal   E[2]={0.0, 0.0};
  PetscReal   vxm, vym, vxt, vxp, vyp;
  PetscInt    ndim=2;

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
  for (dim=0; dim<2; dim++) {
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
    x = pos[ip*ndim + 0] / dx;
    y = pos[ip*ndim + 1] / dy;
    /* Compute the electric field due to this particle: \vec{E} = -\nabla\phi. */
    PetscCall(DifferenceVector2D(phi, x, y, ctx->grid, E));
    vxm = vel[ip*ndim + 0] + Escale[0]*E[0];
    vym = vel[ip*ndim + 1] + Escale[1]*E[1];
    vxt = vxm + vym*t;
    vyp = vym - vxt*s;
    vxp = vxt + vyp*t;
    /* Assign new particle velocities. */
    vel[ip*ndim + 0] = vxp;
    vel[ip*ndim + 1] = vyp;
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
PetscErrorCode BorisMoverBz3D(PetscReal dt, void *opts)
{
  Context    *ctx=(Context *)opts;
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
  PetscInt     ndim=3;

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
    x = pos[ip*ndim + 0] / dx;
    y = pos[ip*ndim + 1] / dy;
    z = pos[ip*ndim + 2] / dz;
    /* Compute the electric field due to this particle: \vec{E} = -\nabla\phi. */
    PetscCall(DifferenceVector3D(phi, x, y, z, ctx->grid, E));
    vxm = vel[ip*ndim + 0] + Escale[0]*E[0];
    vym = vel[ip*ndim + 1] + Escale[1]*E[1];
    vz  = vel[ip*ndim + 2] + Escale[2]*E[2];
    vxt = vxm + vym*t;
    vyp = vym - vxt*s;
    vxp = vxt + vyp*t;
    /* Assign new particle velocities. */
    vel[ip*ndim + 0] = vxp;
    vel[ip*ndim + 1] = vyp;
    vel[ip*ndim + 2] = vz;
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
PetscErrorCode UpdateVelocities(PetscInt ndim, PetscReal dt, Context *ctx)
{
  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Apply a function to integrate dv/dt = E + vxB.
  PetscCall(ctx->ions.push(dt, ctx));

  // Apply the appropriate collision algorithm.
  if (ctx->neutrals.m > 0.0) {
    PetscCall(ctx->ions.collide(ndim, dt, ctx));
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


