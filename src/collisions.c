#include "ecsprt.h"
#include "random.h"

/* Basic 3-D elastic-scattering routine for ion-neutral collisions.

This function is based on, but not identical to, EPPIC elastic_scatter.
*/
PetscErrorCode ScatterElastic(PetscInt ndim, PetscReal dt, void *opts)
{
  Context   *ctx=(Context *)opts;                                              // the application context
  PetscInt   Nc;                                                               // the number of collisions to attempt
  PetscInt   Ns=0;                                                             // the number of successful collision attempts
  PetscInt   Nf=0;                                                             // the number of failed collision attempts
  PetscInt   np;                                                               // the number of ions on this rank
  PetscInt   ip;                                                               // the current ion
  PetscInt   dim;
  PetscReal  fc=ctx->ions.nu * dt;                                             // the product of the collision rate and the time step
  PetscReal  viT=ctx->ions.vT;                                                 // the ion-species thermal speed
  PetscReal  vi0[3]={ctx->ions.v0x, ctx->ions.v0y, ctx->ions.v0z};             // zeroth-order ion velocity
  PetscReal  vn0[3]={ctx->neutrals.v0x, ctx->neutrals.v0y, ctx->neutrals.v0z}; // zeroth-order neutral velocity
  PetscReal  vnT=ctx->neutrals.vT;                                             // the neutral-species thermal speed
  PetscReal  vrm=0.0;                                                          // the maximum ion-neutral relative velocity
  PetscReal  tmp;
  PetscReal  mi=ctx->ions.m;                                                   // the ion-species mass
  PetscReal  mn=ctx->neutrals.m;                                               // the neutral-species mass
  PetscReal  M=mn+mi;                                                          // the total mass (mi+mn)
  PetscReal *vel;
  PetscReal  vn[3];                                                            // neutral-particle velocity
  PetscReal  vi[3];                                                            // ion velocity
  PetscReal  vr[3]={0.0, 0.0, 0.0};                                            // ion-neutral relative-velocity
  PetscReal  vc[3];                                                            // center-of-mass velocity
  PetscReal  vrr;                                                              // ion-neutral relative-velocity magnitude
  PetscReal  vcr;                                                              // the ion speed with respect to the center of mass
  PetscReal  costht, sintht, cosphi, sinphi;
  PetscReal  uperp, uphi;                                                      // components of the unit scattering vector
  PetscReal  u0x, u0y, u0z;
  PetscReal  u1x, u1y, u1z;
  PetscReal  u2[3], vf[3], vfr;
  long       seed=getseed(*ctx);
  PetscReal  ratio;                                                            // ratio of current ion's final speed to thermal speed

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Compute the maximum relative velocity.
  tmp = 0.0;
  for (dim=0; dim<ndim; dim++) {
    tmp += PetscSqr(vi0[dim] - vn0[dim]);
  }
  vrm = 4.0*viT + PetscSqrtReal(tmp);

  // Get the number of ions on this rank.
  PetscCall(DMSwarmGetLocalSize(ctx->swarmDM, &np));

  // Compute the number of collisions to attempt.
  Nc = (PetscInt)(((PetscReal)np * fc * vrm) / ((vnT + ctx->neutrals.v0) * PetscSqrtReal(mn / mi)));
  if (Nc > np) {
    Nc = np;
  } else if (Nc < 0) {
    Nc = 0;
  }
  ctx->log.ranks("[%d] Colliding %d particles out of %d ...\n", ctx->mpi.rank, Nc, np);

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(ctx->swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Attempt collisions until we reach the required number.
  while (Ns < Nc) {
    // Choose a random ion from the full distribution.
    ip = (PetscInt)(np*ran3(&seed));

    // - Store the ion velocity components.
    // - Choose a random neutral-particle velocity.
    // - Compute the ion-neutral relative velocity components and magnitude.
    tmp = 0.0;
    for (dim=0; dim<ndim; dim++) {
      vi[dim] = vel[ip*ndim + dim];
      vn[dim] = gasdev(&seed)*vnT + vn0[dim];
      vr[dim] = vi[dim] - vn[dim];
      tmp += PetscSqr(vr[dim]);
    }
    vrr = PetscSqrtReal(tmp);

    // Collide only if the squared relative velocity is greater than the square
    // of a random percentage of the maximum relative velocity.
    if (PetscSqr(ran3(&seed)) < PetscSqr(vrr) / PetscSqr(vrm)) {

      // - Compute the center-of-mass (CoM) velocity.
      // - Compute the particle speed relative to the CoM velocity.
      tmp = 0.0;
      for (dim=0; dim<ndim; dim++) {
        vc[dim] = (vi[dim]*mi + vn[dim]*mn) / M;
        tmp += PetscSqr(vi[dim] - vc[dim]);
      }
      vcr = PetscSqrtReal(tmp);

      // Compute the unit scattering vector relative to the CoM z axis.
      u0z = 2.0*ran3(&seed) - 1.0;
      uperp = PetscSqrtReal(1.0 - PetscSqr(u0z));
      uphi = 2*PETSC_PI * ran3(&seed);
      u0x = uperp*PetscCosReal(uphi);
      u0y = uperp*PetscSinReal(uphi);

      // Rotate the CoM frame to align its z axis with the incident direction.
      costht = vr[2] / vrr;
      sintht = PetscSqrtReal(1.0 - PetscSqr(costht));
      cosphi = vr[0] / (vrr*sintht);
      sinphi = vr[1] / (vrr*sintht);

      /* Rotate the unit scattering vector to the incident coordinate system.
      1. rotation about CoM y axis:          (xc, yc, zc) -> (xp, yp, zp)
      2. rotation about intermediate z axis: (xp, yp, zp) -> (xi, yi, zi)
      */
      u1x = u0z*sintht + u0x*costht;
      u1y = u0y;
      u1z = u0z*costht - u0x*sintht;
      u2[0] = u1x*cosphi - u1y*sinphi;
      u2[1] = u1x*sinphi + u1y*cosphi;
      u2[2] = u1z;

      /* Assign final CoM velocity components.
      vfx = vcr * ((uz*sintht + ux*costht)*cosphi - uy*sinphi)
      vfy = vcr * ((uz*sintht - ux*costht)*sinphi + uy*cosphi)
      vfz = vcr * ( uz*costht - ux*sintht                    )
      */
      tmp = 0.0;
      for (dim=0; dim<ndim; dim++) {
        vf[dim] = vc[dim] + vcr*u2[dim];
        tmp += PetscSqr(vf[dim]);
      }

      /* Finalize
      - if result is unphysical, do not use it
      - otherwise, update the ion velocity components and count the collision
      */
      vfr = PetscSqrtReal(tmp);
      ratio = vfr / viT;
      ctx->log.ranks("vfr / viT = %6.4f / %6.4f = %6.4f\n", vfr, viT, ratio);
      if (ratio > 10) {
        ctx->log.self("[%d] Warning: Refusing to accept collision that results in final speed = %4.1f times thermal speed\n", ctx->mpi.rank, ratio);
        Nf++;
        // Terminate the simulation if at least 10 collisions have failed.
        if (Nf >= 10) {
          ctx->log.status("Failed to collide %d ion-neutral pairs. Aborting.\n\n", Nf);
          MPI_Abort(PETSC_COMM_WORLD, 1);
        }
      } else {
        for (dim=0; dim<ndim; dim++) {
          vel[ip*ndim + dim] = vf[dim];
        }
        Ns++;
      }

    }
  }

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(ctx->swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Report the number of actual collisions.
  ctx->log.world("Collision efficiency: %6.4f.\n", (PetscReal)Ns/(PetscReal)(Ns+Nf));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


