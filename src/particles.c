#include <petsc.h>
#include "ecsprt.h"
#include "positions.h"
#include "velocities.h"
#include "random.h"
#include "vector-math.h"
#include "particles.h"
#include "pic.h"

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


PetscErrorCode Apply2DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  Context   *ctx=(Context *)opts;
  PetscInt   ip, ndim=2;
  PetscReal  x, y;
  PetscReal  Lx=ctx->grid.Lx;
  PetscReal  Ly=ctx->grid.Ly;
  PetscReal  x0=ctx->grid.x0;
  PetscReal  y0=ctx->grid.y0;
  PetscReal  x1=ctx->grid.x1;
  PetscReal  y1=ctx->grid.y1;

  PetscFunctionBeginUser;

  // TODO: Implement other BC in ion loop. Some of those BC will require
  // modifying the velocity (e.g., vx = -vx for reflection at the x boundary, or
  // vx = vx0 for injection/advection along the x axis).

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    // Extract current positions.
    x = pos[ip*ndim + 0];
    y = pos[ip*ndim + 1];
    // Update the x position.
    if (ctx->ions.xBC == BC_PERIODIC) {
      if (x < x0) {x += Lx;}
      if (x >= x1) {x -= Lx;}
    }
    // Update the y position.
    if (ctx->ions.yBC == BC_PERIODIC) {
      if (y < y0) {y += Ly;}
      if (y >= y1) {y -= Ly;}
    }
    // Copy new positions.
    pos[ip*ndim + 0] = x;
    pos[ip*ndim + 1] = y;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode Apply3DBC(PetscInt np, PetscReal *pos, PetscReal *vel, void *opts)
{
  Context   *ctx=(Context *)opts;
  PetscInt  ip, ndim=3;
  PetscReal x, y, z;
  PetscReal Lx=ctx->grid.Lx;
  PetscReal Ly=ctx->grid.Ly;
  PetscReal Lz=ctx->grid.Lz;
  PetscReal x0=ctx->grid.x0;
  PetscReal y0=ctx->grid.y0;
  PetscReal z0=ctx->grid.z0;
  PetscReal x1=ctx->grid.x1;
  PetscReal y1=ctx->grid.y1;
  PetscReal z1=ctx->grid.z1;

  PetscFunctionBeginUser;

  // TODO: Implement other BC in ion loop. Some of those BC will require
  // modifying the velocity (e.g., vx = -vx for reflection at the x boundary, or
  // vx = vx0 for injection/advection along the x axis).

  // Loop over ions.
  for (ip=0; ip<np; ip++) {
    // Extract current positions.
    x = pos[ip*ndim + 0];
    y = pos[ip*ndim + 1];
    z = pos[ip*ndim + 2];
    // Update the x position.
    if (ctx->ions.xBC == BC_PERIODIC) {
      if (x < x0) {x += Lx;}
      if (x >= x1) {x -= Lx;}
    }
    // Update the y position.
    if (ctx->ions.yBC == BC_PERIODIC) {
      if (y < y0) {y += Ly;}
      if (y >= y1) {y -= Ly;}
    }
    // Update the z position.
    if (ctx->ions.zBC == BC_PERIODIC) {
      if (z < z0) {z += Lz;}
      if (z >= z1) {z -= Lz;}
    }
    // Copy new positions.
    pos[ip*ndim + 0] = x;
    pos[ip*ndim + 1] = y;
    pos[ip*ndim + 2] = z;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Enforce boundary conditions and migrate particles.

This routine should be the only place where the code calls `DMSwarmMigrate`.
*/
PetscErrorCode ApplyBCAndMigrate(Context *ctx)
{
  DM         swarmDM=ctx->swarmDM;
  PetscReal *pos, *vel;
  PetscInt   np, Np0, Np1;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Get an array representation of the ion positions.
  PetscCall(DMSwarmGetField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Get an array representation of the ion velocities.
  PetscCall(DMSwarmGetField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Get the number of particles on this rank.
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));

  // Apply boundary conditions.
  PetscCall(ctx->ions.applyBC(np, pos, vel, ctx));

  // Restore the ion-velocities array.
  PetscCall(DMSwarmRestoreField(swarmDM, "velocity", NULL, NULL, (void **)&vel));

  // Restore the ion-positions array.
  PetscCall(DMSwarmRestoreField(swarmDM, DMSwarmPICField_coor, NULL, NULL, (void **)&pos));

  // Echo pre-migration sizes.
  PetscCall(DMSwarmGetSize(swarmDM, &Np0));
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));
  ctx->log.world("\n");
  ctx->log.ranks("[%d] Local # of ions before migration: %d\n", ctx->mpi.rank, np);
  ctx->log.world("   Global # of ions before migration: %d\n", Np0);

  // Migrate ions among processes.
  PetscCall(DMSwarmMigrate(swarmDM, PETSC_TRUE));

  // Echo post-migration sizes.
  PetscCall(DMSwarmGetSize(swarmDM, &Np1));
  PetscCall(DMSwarmGetLocalSize(swarmDM, &np));
  ctx->log.world("\n");
  ctx->log.ranks("[%d] Local # of ions after migration: %d\n", ctx->mpi.rank, np);
  ctx->log.world("   Global # of ions after migration: %d\n", Np1);

  if ((ctx->ions.xBC == BC_PERIODIC) && (ctx->ions.yBC == BC_PERIODIC) && (ctx->ions.zBC == BC_PERIODIC)) {
    if (Np1 != Np0) {
      ctx->log.status("\n");
      ctx->log.status("ERROR: Total number of ions has changed (%+d).\n", Np1-Np0);
      ctx->log.status("       This should not happen with fully periodic boundary conditions.\n");
      ctx->log.status("\n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


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


/* Compute moments of the ion distribution. */
PetscErrorCode CollectFluidMoments(Context *ctx)
{
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
    x = pos[ip*NDIM + 0] / dx;
    y = pos[ip*NDIM + 1] / dy;
    z = pos[ip*NDIM + 2] / dz;

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
    for (dim=0, dof=1; dim<NDIM; dim++, dof++) {
      w = vel[ip*NDIM + dim];
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


/* Apply the standard 3-D Boris-mover algorithm to ion velocities.

See latter part of Birdsall & Langdon section 4-4.
*/
PetscErrorCode BorisMover3D(PetscReal dt, Context *ctx)
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


/* Apply the standard 1-D Boris-mover algorithm to ion velocities.

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


