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


/* Basic elastic-scattering routine for ion-neutral collisions.

This function is based on, but not identical to, EPPIC elastic_scatter.
*/
PetscErrorCode ComputeCollisions(PetscReal dt, Context *ctx)
{
  PetscInt   Nc;                                      // the number of collisions to attempt
  PetscInt   Ns=0;                                    // the number of successful collision attempts
  PetscInt   Nf=0;                                    // the number of failed collision attempts
  PetscInt   np;                                      // the number of ions on this rank
  PetscInt   ip;                                      // the current ion
  PetscReal  fc=ctx->ions.nu * dt;                    // the product of the collision rate and the time step
  PetscReal  viT=ctx->ions.vT;                        // the ion-species thermal speed
  PetscReal  vi0x=ctx->ions.v0x;                      // the ion-species x-axis drift component
  PetscReal  vi0y=ctx->ions.v0y;                      // the ion-species y-axis drift component
  PetscReal  vi0z=ctx->ions.v0z;                      // the ion-species z-axis drift component
  PetscReal  vnT=ctx->neutrals.vT;                    // the neutral-species thermal speed
  PetscReal  vn0=ctx->neutrals.v0;                    // the neutral-species drift speed
  PetscReal  vn0x=ctx->neutrals.v0x;                  // the neutral-species x-axis drift component
  PetscReal  vn0y=ctx->neutrals.v0y;                  // the neutral-species y-axis drift component
  PetscReal  vn0z=ctx->neutrals.v0z;                  // the neutral-species z-axis drift component
  PetscReal  vrm;                                     // the maximum ion-neutral relative velocity
  PetscReal  mi=ctx->ions.m;                          // the ion-species mass
  PetscReal  mn=ctx->neutrals.m;                      // the neutral-species mass
  PetscReal  M=mn+mi;                                 // the total mass (mi+mn)
  DM         swarmDM=ctx->swarmDM;
  PetscReal *vel;
  PetscReal  vnx, vny, vnz;                           // neutral-particle velocity components
  PetscReal  vix, viy, viz;                           // ion velocity components
  PetscReal  vrx, vry, vrz;                           // ion-neutral relative-velocity components
  PetscReal  vrr;                                     // ion-neutral relative-velocity magnitude
  PetscReal  vcx, vcy, vcz;                           // center-of-mass velocity components
  PetscReal  vcr;                                     // the ion speed with respect to the center of mass
  PetscReal  costht, sintht, cosphi, sinphi;
  PetscReal  ux, uy, uz, uperp, uphi;                 // components of the unit scattering vector
  PetscReal  ux1, uy1, uz1, ux2, uy2, uz2;
  PetscReal  vfx, vfy, vfz, vfr;                      // components and magnitude of current ion's final speed
  long       seed=getseed(*ctx);
  PetscReal  ratio;                                   // ratio of current ion's final speed to thermal speed

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Compute the maximum relative velocity.
  vrm = 4.0*viT + PetscSqrtReal(PetscSqr(vi0x-vn0x) + PetscSqr(vi0y-vn0y) + PetscSqr(vi0z-vn0z));

  // Get the number of ions on this rank.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Compute the number of collisions to attempt.
  Nc = (PetscInt)(((PetscReal)np * fc * vrm) / ((vnT + vn0) * PetscSqrtReal(mn / mi)));
  if (Nc > np) {
    Nc = np;
  } else if (Nc < 0) {
    Nc = 0;
  }
  ctx->log.ranks("[%d] Colliding %d particles out of %d ...\n", ctx->mpi.rank, Nc, np);

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Attempt collisions until we reach the required number.
  while (Ns < Nc) {
    // Choose a random ion from the full distribution.
    ip = (PetscInt)(np*ran3(&seed));

    // Store the ion velocity components.
    vix = vel[ip*NDIM + 0];
    viy = vel[ip*NDIM + 1];
    viz = vel[ip*NDIM + 2];

    // Choose a random neutral-particle velocity.
    vnx = gasdev(&seed)*vnT + vn0x;
    vny = gasdev(&seed)*vnT + vn0y;
    vnz = gasdev(&seed)*vnT + vn0z;

    // Compute the ion-neutral relative velocity.
    vrx = vix - vnx;
    vry = viy - vny;
    vrz = viz - vnz;
    vrr = PetscSqrtReal(PetscSqr(vrx) + PetscSqr(vry) + PetscSqr(vrz));

    // Collide only if the squared relative velocity is greater than the square
    // of a random percentage of the maximum relative velocity.
    if (PetscSqr(ran3(&seed)) < PetscSqr(vrr) / PetscSqr(vrm)) {

      // Compute the center-of-mass (CoM) velocity.
      vcx = (vix*mi + vnx*mn) / M;
      vcy = (viy*mi + vny*mn) / M;
      vcz = (viz*mi + vnz*mn) / M;

      // Compute the particle speed relative to the CoM velocity.
      vcr = PetscSqrtReal(PetscSqr(vix-vcx) + PetscSqr(viy-vcy) + PetscSqr(viz-vcz));

      // Compute the unit scattering vector relative to the CoM z axis.
      uz = 2.0*ran3(&seed) - 1.0;
      uperp = PetscSqrtReal(1.0 - PetscSqr(uz));
      uphi = 2*PETSC_PI * ran3(&seed);
      ux = uperp*PetscCosReal(uphi);
      uy = uperp*PetscSinReal(uphi);

      // Rotate the CoM frame to align its z axis with the incident direction.
      costht = vrz / vrr;
      sintht = PetscSqrtReal(1.0 - PetscSqr(costht));
      cosphi = vrx / (vrr*sintht);
      sinphi = vry / (vrr*sintht);

      /* Rotate the unit scattering vector to the incident coordinate system.
      1. rotation about CoM y axis:          (xc, yc, zc) -> (xp, yp, zp)
      2. rotation about intermediate z axis: (xp, yp, zp) -> (xi, yi, zi)
      */
      ux1 = uz*sintht + ux*costht;
      uy1 = uy;
      uz1 = uz*costht - ux*sintht;
      ux2 = ux1*cosphi - uy1*sinphi;
      uy2 = ux1*sinphi + uy1*cosphi;
      uz2 = uz1;

      /* Assign final CoM velocity components.
      vfx = vcr * ((uz*sintht + ux*costht)*cosphi - uy*sinphi)
      vfy = vcr * ((uz*sintht - ux*costht)*sinphi + uy*cosphi)
      vfz = vcr * ( uz*costht - ux*sintht                    )
      */
      vfx = vcx + vcr*ux2;
      vfy = vcy + vcr*uy2;
      vfz = vcz + vcr*uz2;

      /* Finalize
      - if result is unphysical, do not use it
      - otherwise, update the ion velocity components and count the collision
      */
      vfr = PetscSqrtReal(PetscSqr(vfx) + PetscSqr(vfy) + PetscSqr(vfz));
      ratio = vfr / viT;
      if (ratio > 10) {
        ctx->log.self("[%d] Warning: Refusing to accept collision that results in final speed = %4.1f times thermal speed\n", ctx->mpi.rank, ratio);
        Nf++;
        // Terminate the simulation if at least 10 collisions have failed.
        if (Nf >= 10) {
          ctx->log.status("Failed to collide %d ion-neutral pairs. Aborting.\n\n", Nf);
          MPI_Abort(PETSC_COMM_WORLD, 1);
        }
      } else {
        vel[ip*NDIM + 0] = vfx;
        vel[ip*NDIM + 1] = vfy;
        vel[ip*NDIM + 2] = vfz;
        Ns++;
      }

    }
  }

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Report the number of actual collisions.
  ctx->log.world("Collision efficiency: %6.4f.\n", (PetscReal)Ns/(PetscReal)(Ns+Nf));

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
    PetscCall(ComputeCollisions(dt, ctx));
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


