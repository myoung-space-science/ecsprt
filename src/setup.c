#include <petsc.h>
#include "ecsprt.h"
#include "lhs.h"
#include "rhs.h"
#include "logging.h"
#include "moments.h"
#include "boundaries.h"
#include "velocities.h"
#include "collisions.h"


/* Set parameter values common to simulation and solver applications. */
PetscErrorCode SetUpContext(CLI cli, Context *ctx)
{

  PetscReal tmp;

  PetscFunctionBeginUser;

  // Declare the name of the options log.
  PetscCall(PetscStrcpy(ctx->optionsLog, "options.log"));

  // Set the LHS function based on LHS type.
  switch (cli.lhsType) {
  case LHS_IDENTITY:
    ctx->potential.lhs = ComputeIdentityLHS;
    ctx->potential.stencilSize = 1;
    ctx->potential.stencilType = DMDA_STENCIL_STAR;
    break;
  case LHS_LAPLACIAN:
    switch (cli.ndim) {
    case 2:
      ctx->potential.lhs = ComputeLaplacianLHS2D;
      ctx->potential.stencilSize = 5;
      break;
    case 3:
      ctx->potential.lhs = ComputeLaplacianLHS3D;
      ctx->potential.stencilSize = 7;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", cli.ndim);
    }
    ctx->potential.stencilType = DMDA_STENCIL_STAR;
    break;
  case LHS_FULL:
    switch (cli.ndim) {
    case 2:
      ctx->potential.lhs = ComputeFullLHS2D;
      ctx->potential.stencilSize = 9;
      break;
    case 3:
      ctx->potential.lhs = ComputeFullLHS3D;
      ctx->potential.stencilSize = 11;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", cli.ndim);
    }
    ctx->potential.stencilType = DMDA_STENCIL_BOX;
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown LHS type: \"%s\"\n", LHSTypes[cli.lhsType]);
  }

  // Set RHS function based on RHS type.
  switch (cli.rhsType) {
  case RHS_CONSTANT:
    ctx->potential.rhs = ComputeConstantRHS;
    break;
  case RHS_SINUSOIDAL:
    switch (cli.ndim)
    {
    case 2:
      ctx->potential.rhs = ComputeSinusoidalRHS2D;
      break;
    case 3:
      ctx->potential.rhs = ComputeSinusoidalRHS3D;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", cli.ndim);
    }
    break;
  case RHS_FULL:
    switch (cli.ndim)
    {
    case 2:
      ctx->potential.rhs = ComputeFullRHS2D;
      break;
    case 3:
      ctx->potential.rhs = ComputeFullRHS3D;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", cli.ndim);
    }
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown RHS type: \"%s\"\n", RHSTypes[cli.rhsType]);
  }

  switch (cli.ndim)
  {
  case 2:
    ctx->ions.applyBC = Apply2DBC;
    ctx->ions.collect = Collect2DFluidMoments;
    ctx->ions.push    = BorisMoverBz2D;
    break;
  case 3:
    ctx->ions.applyBC = Apply3DBC;
    ctx->ions.collect = Collect3DFluidMoments;
    ctx->ions.push    = BorisMoverBz3D;
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", cli.ndim);
  }
  ctx->ions.collide = ScatterElastic;

  // Set fundamental parameter values.
  ctx->electrons.q = -Q;
  ctx->electrons.m = ME;

  // Explicitly set the neutral charge to zero.
  ctx->neutrals.q = 0.0;

  // Copy plasma and species parameter values.
  ctx->plasma.n0     = cli.n0;
  ctx->plasma.B0     = cli.B0;
  ctx->plasma.E0     = cli.E0;
  ctx->plasma.Np     = cli.Np;
  ctx->ions.q        = cli.qi;
  ctx->ions.m        = cli.mi;
  ctx->ions.v0x      = cli.vi0x;
  ctx->ions.v0y      = cli.vi0y;
  ctx->ions.v0z      = cli.vi0z;
  ctx->ions.v0       = cli.vi0;
  ctx->ions.vTx      = cli.viTx;
  ctx->ions.vTy      = cli.viTy;
  ctx->ions.vTz      = cli.viTz;
  ctx->ions.vT       = cli.viT;
  ctx->ions.T        = cli.Ti;
  ctx->ions.nu       = cli.nui;
  ctx->electrons.v0x = cli.ve0x;
  ctx->electrons.v0y = cli.ve0y;
  ctx->electrons.v0z = cli.ve0z;
  ctx->electrons.v0  = cli.ve0;
  ctx->electrons.vTx = cli.veTx;
  ctx->electrons.vTy = cli.veTy;
  ctx->electrons.vTz = cli.veTz;
  ctx->electrons.vT  = cli.veT;
  ctx->electrons.T   = cli.Te;
  ctx->electrons.nu  = cli.nue;
  ctx->neutrals.m    = cli.mn;
  ctx->neutrals.v0x  = cli.vn0x;
  ctx->neutrals.v0y  = cli.vn0y;
  ctx->neutrals.v0z  = cli.vn0z;
  ctx->neutrals.v0   = cli.vn0;
  ctx->neutrals.vTx  = cli.vnTx;
  ctx->neutrals.vTy  = cli.vnTy;
  ctx->neutrals.vTz  = cli.vnTz;
  ctx->neutrals.vT   = cli.vnT;
  ctx->neutrals.T    = cli.Tn;

  /* Copy grid parameters. */
  // x dimension
  ctx->grid.Nx = cli.Nx;
  ctx->grid.dx = cli.dx;
  if (cli.x1 == cli.x0) {
      ctx->log.world("Warning: zero-width x dimension\n");
  }
  ctx->grid.x0 = cli.x0;
  ctx->grid.x1 = cli.x1;
  // y dimension
  ctx->grid.Ny = cli.Ny;
  ctx->grid.dy = cli.dy;
  if (cli.y1 == cli.y0) {
      ctx->log.world("Warning: zero-width y dimension\n");
  }
  ctx->grid.y0 = cli.y0;
  ctx->grid.y1 = cli.y1;
  // z dimension
  if (cli.ndim == 3) {
    ctx->grid.Nz = cli.Nz;
    ctx->grid.dz = cli.dz;
    if (cli.z1 == cli.z0) {
        ctx->log.world("Warning: zero-width z dimension\n");
    }
    ctx->grid.z0 = cli.z0;
    ctx->grid.z1 = cli.z1;
  } else {
    ctx->grid.Nz = 1;
    ctx->grid.dz = 0.0;
  }
  /* Set up boundary conditions.
  - If one boundary type for a given dimension is periodic, the other must be
    periodic. Otherwise, the notion of periodicity is meaningless.
  - The fluid (`grid`) boundary condition is periodic if both boundary types are
    periodic, and ghosted if they are not.
  - The specific pair of boundary types for a given dimension, assuming they
    pass the periodicity check, map to a pre-defined particle (`swarm`) boundary
    condition for that dimension.
  */
  // x dimension
  ctx->ions.xBC = cli.xBC;
  if (ctx->ions.xBC == BC_PERIODIC) {
    ctx->grid.xBC = DM_BOUNDARY_PERIODIC;
  } else {
    ctx->grid.xBC = DM_BOUNDARY_GHOSTED;
  }
  // y dimension
  ctx->ions.yBC = cli.yBC;
  if (ctx->ions.yBC == BC_PERIODIC) {
    ctx->grid.yBC = DM_BOUNDARY_PERIODIC;
  } else {
    ctx->grid.yBC = DM_BOUNDARY_GHOSTED;
  }
  // z dimension
  ctx->ions.zBC = cli.zBC;
  if (ctx->ions.zBC == BC_PERIODIC) {
    ctx->grid.zBC = DM_BOUNDARY_PERIODIC;
  } else {
    ctx->grid.zBC = DM_BOUNDARY_GHOSTED;
  }

  // Set species gyrofrequency from q, B0, and m.
  ctx->electrons.Omega = PetscAbsReal(ctx->electrons.q * ctx->plasma.B0 / ctx->electrons.m);
  ctx->ions.Omega = PetscAbsReal(ctx->ions.q * ctx->plasma.B0 / ctx->ions.m);

  // Set species magnetization from Omega and nu.
  ctx->electrons.kappa = ctx->electrons.Omega / ctx->electrons.nu;
  ctx->ions.kappa = ctx->ions.Omega / ctx->ions.nu;

  // Declare scale factor for potential equation.
  ctx->potential.scale = 1.0 / (1 + ctx->electrons.kappa*ctx->electrons.kappa);

  // Compute drift-velocity magnitudes.
  ctx->electrons.v0 = PetscSqrtReal(PetscSqr(ctx->electrons.v0x) + PetscSqr(ctx->electrons.v0y) + PetscSqr(ctx->electrons.v0z));
  ctx->ions.v0 = PetscSqrtReal(PetscSqr(ctx->ions.v0x) + PetscSqr(ctx->ions.v0y) + PetscSqr(ctx->ions.v0z));
  ctx->neutrals.v0 = PetscSqrtReal(PetscSqr(ctx->neutrals.v0x) + PetscSqr(ctx->neutrals.v0y) + PetscSqr(ctx->neutrals.v0z));

  // Make electron temperature and thermal velocity consistent.
  if ((ctx->electrons.vTx != 0.0) || (ctx->electrons.vTy != 0.0) || (ctx->electrons.vTz != 0.0)) {
    ctx->electrons.vT = PetscSqrtReal(PetscSqr(ctx->electrons.vTx) + PetscSqr(ctx->electrons.vTy) + PetscSqr(ctx->electrons.vTz));
    ctx->electrons.T = (0.5 * ctx->electrons.m / KB) * (PetscSqr(ctx->electrons.vT));
  } else {
    tmp = PetscSqrtReal(2.0 * KB * ctx->electrons.T / ctx->electrons.m) / 3.0;
    ctx->electrons.vTx = tmp;
    ctx->electrons.vTy = tmp;
    ctx->electrons.vTz = tmp;
    ctx->electrons.vT = PetscSqrtReal(PetscSqr(ctx->electrons.vTx) + PetscSqr(ctx->electrons.vTy) + PetscSqr(ctx->electrons.vTz));
  }

  // Make ion temperature and thermal velocity consistent.
  if ((ctx->ions.vTx != 0.0) || (ctx->ions.vTy != 0.0) || (ctx->ions.vTz != 0.0)) {
    ctx->ions.vT = PetscSqrtReal(PetscSqr(ctx->ions.vTx) + PetscSqr(ctx->ions.vTy) + PetscSqr(ctx->ions.vTz));
    ctx->ions.T = (0.5 * ctx->ions.m / KB) * (PetscSqr(ctx->ions.vT));
  } else {
    tmp = PetscSqrtReal(2.0 * KB * ctx->ions.T / ctx->ions.m) / 3.0;
    ctx->ions.vTx = tmp;
    ctx->ions.vTy = tmp;
    ctx->ions.vTz = tmp;
    ctx->ions.vT = PetscSqrtReal(PetscSqr(ctx->ions.vTx) + PetscSqr(ctx->ions.vTy) + PetscSqr(ctx->ions.vTz));
  }

  // Set default neutral temperature equal to ion temperature.
  if (ctx->neutrals.T == -1.0) {
    ctx->neutrals.T = ctx->ions.T;
    ctx->log.world("Warning: Setting neutral temperature equal to ion temperature (%.1f K)\n", ctx->ions.T);
  }

  // Set neutral thermal velocity from temperature.
  ctx->neutrals.vT = PetscSqrtReal(2.0 * KB * ctx->neutrals.T / ctx->neutrals.m);

  // [DEV] Hard-code the electron thermal coefficient to 1.0 (isothermal).
  ctx->gammaT = 1.0;

  // TODO: Should we set default collision frequencies based on an analytic
  // formulation (e.g., from Schunk & Nagy)?

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Destroy the application context and associated objects. */
PetscErrorCode DestroyContext(Context *ctx)
{

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(VecDestroy(&ctx->moments));
  PetscCall(DMDestroy(&ctx->fluidDM));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode NormalizeGrid(DM dm, Context *ctx)
{
  PetscInt Nx, Ny, Nz;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Synchronize values of Nx, Ny, and Nz.
  PetscCall(DMDAGetInfo(dm, NULL, &Nx, &Ny, &Nz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
  if (ctx->grid.Nx == -1) {
    ctx->grid.Nx = Nx;
  }
  if (ctx->grid.Ny == -1) {
    ctx->grid.Ny = Ny;
  }
  if (ctx->grid.Nz == -1) {
    ctx->grid.Nz = Nz;
  }
  /* Set the physical grid attributes.

  At this point, the code has set the values of Nx, Ny, and Nz. Let q represent
  each of {x,y,z}. If the user did not pass in a positive value for dq, this
  section will set Lq = q1-q0 and compute dq in terms of Lq. Otherwise it will
  compute Lq, q0, and q1 in terms of dq and Nq. Note that the latter case may
  result in overwriting user values of q0 and q1.
  */
  if (ctx->grid.dx <= 0.0) {
    ctx->grid.Lx = ctx->grid.x1 - ctx->grid.x0;
    ctx->grid.dx = ctx->grid.Lx / (PetscReal)ctx->grid.Nx;
  } else {
    ctx->grid.Lx = ctx->grid.dx * (PetscReal)ctx->grid.Nx;
    ctx->grid.x0 = 0.0;
    ctx->grid.x1 = ctx->grid.Lx;
  }
  if (ctx->grid.dy <= 0.0) {
    ctx->grid.Ly = ctx->grid.y1 - ctx->grid.y0;
    ctx->grid.dy = ctx->grid.Ly / (PetscReal)ctx->grid.Ny;
  } else {
    ctx->grid.Ly = ctx->grid.dy * (PetscReal)ctx->grid.Ny;
    ctx->grid.y0 = 0.0;
    ctx->grid.y1 = ctx->grid.Ly;
  }
  if (ctx->grid.dz <= 0.0) {
    ctx->grid.Lz = ctx->grid.z1 - ctx->grid.z0;
    ctx->grid.dz = ctx->grid.Lz / (PetscReal)ctx->grid.Nz;
  } else {
    ctx->grid.Lz = ctx->grid.dz * (PetscReal)ctx->grid.Nz;
    ctx->grid.z0 = 0.0;
    ctx->grid.z1 = ctx->grid.Lz;
  }

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Create the data manager for Eulerian fluid quantities.

This routine ultimately determines the number of grid cells in each dimension
(Nx, Ny, and Nz) and the total number of charged particles (Np). In the former
case, we want to let the user specify Nx, Ny, or Nz via the PETSc options
`-da_grid_x`, `-da_grid_y`, and `-da_grid_z`. In the latter case, we want the
default value of Np to correspond to one particle per cell, which we can't
compute until we are certain about the values of Nx, Ny, and Nz.
*/
PetscErrorCode CreateFluidDM(PetscInt ndim, Context *ctx)
{
  PetscInt        Nx=(ctx->grid.Nx > 0 ? ctx->grid.Nx : 7);
  PetscInt        Ny=(ctx->grid.Ny > 0 ? ctx->grid.Ny : 7);
  PetscInt        Nz=(ctx->grid.Nz > 0 ? ctx->grid.Nz : 7);
  DMBoundaryType  xBC=ctx->grid.xBC;
  DMBoundaryType  yBC=ctx->grid.yBC;
  DMBoundaryType  zBC=ctx->grid.zBC;
  DMDAStencilType stencilType=ctx->potential.stencilType;
  PetscInt        dof=4;
  PetscInt        width=1;
  DM              dm;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Create the DM object in 2 or 3 dimensions.
  switch (ndim)
  {
  case 2:
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, xBC, yBC, stencilType, Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, &dm));
    break;
  case 3:
    PetscCall(DMDACreate3d(PETSC_COMM_WORLD, xBC, yBC, zBC, stencilType, Nx, Ny, Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, NULL, &dm));
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", ndim);
  }
  // Perform basic setup.
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)dm, "fluid_"));
  PetscCall(DMDASetElementType(dm, DMDA_ELEMENT_Q1));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(PetscObjectSetName((PetscObject)dm, "FluidDM"));
  // Synchronize values of grid parameters.
  PetscCall(NormalizeGrid(dm, ctx));
  // Set uniform coordinates on the DM.
  PetscCall(DMDASetUniformCoordinates(dm, ctx->grid.x0, ctx->grid.x1, ctx->grid.y0, ctx->grid.y1, ctx->grid.z0, ctx->grid.z1));
  // Declare grid-quantity names.
  PetscCall(DMDASetFieldName(dm, 0, "density"));
  PetscCall(DMDASetFieldName(dm, 1, "x flux"));
  PetscCall(DMDASetFieldName(dm, 2, "y flux"));
  PetscCall(DMDASetFieldName(dm, 3, "z flux"));
  // Create a persistent vector for outputing fluid quantities.
  PetscCall(DMCreateGlobalVector(dm, &ctx->moments));
  PetscCall(VecZeroEntries(ctx->moments));
  // Assign the grid DM to the application context.
  ctx->fluidDM = dm;

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Create the data manager for the electrostatic potential.

The logical grid for the electrostatic potential will have the same shape,
stencil type, and boundary conditions as the logical grid for fluid quantities.
Unlike the latter, it will use a stencil width of 2 in order to accommodate
second-order finite differencing of the electric field.

This function should follow `CreateFluidDM`, which performs various set-up and
sychronization tasks on grid parameters.
*/
PetscErrorCode CreatePotentialDM(PetscInt ndim, Context *ctx)
{
  PetscInt        Nx=ctx->grid.Nx;
  PetscInt        Ny=ctx->grid.Ny;
  PetscInt        Nz=ctx->grid.Nz;
  DMBoundaryType  xBC=ctx->grid.xBC;
  DMBoundaryType  yBC=ctx->grid.yBC;
  DMBoundaryType  zBC=ctx->grid.zBC;
  DMDAStencilType stencilType=ctx->potential.stencilType;
  PetscInt        dof=1;
  PetscInt        width=2;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Create the DM object in 2 or 3 dimensions.
  switch (ndim)
  {
  case 2:
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, xBC, yBC, stencilType, Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, &ctx->potential.dm));
    break;
  case 3:
    PetscCall(DMDACreate3d(PETSC_COMM_WORLD, xBC, yBC, zBC, stencilType, Nx, Ny, Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, NULL, &ctx->potential.dm));
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unsupported spatial dimension: %d", ndim);
  }
  // Perform basic setup.
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)(ctx->potential.dm), "potential_"));
  PetscCall(DMDASetElementType(ctx->potential.dm, DMDA_ELEMENT_Q1));
  PetscCall(DMSetFromOptions(ctx->potential.dm));
  PetscCall(DMSetUp(ctx->potential.dm));
  PetscCall(PetscObjectSetName((PetscObject)(ctx->potential.dm), "PotentialDM"));
  // Set uniform coordinates on the DM.
  PetscCall(DMDASetUniformCoordinates(ctx->potential.dm, ctx->grid.x0, ctx->grid.x1, ctx->grid.y0, ctx->grid.y1, ctx->grid.z0, ctx->grid.z1));
  // Assign the field name.
  PetscCall(DMDASetFieldName(ctx->potential.dm, 0, "potential"));
  // Associate the user context with this DM.
  PetscCall(DMSetApplicationContext(ctx->potential.dm, &ctx));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Create the data manager for ion distributions.

The logical grid for the ion distributions will have the same shape, stencil
width, degrees of freedom, and boundary conditions as the logical grid for fluid
quantities. Unlike the latter, it will unconditionally use a box stencil type in
order to accommodate nearest-neighbor moment collection.

This function calls `CreateFluidDM`, which performs various set-up and
sychronization tasks on grid parameters.
*/
PetscErrorCode CreateSwarmDM(PetscInt ndim, Context *ctx)
{
  PetscInt        dim;
  DM              swarmDM;
  PetscInt        bufsize=0;
  PetscInt        np;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Create the fluid DM in 2 or 3 dimensions.
  PetscCall(CreateFluidDM(ndim, ctx));

  // Create the swarm DM.
  PetscCall(DMCreate(PETSC_COMM_WORLD, &swarmDM));
  // Perform basic setup.
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)swarmDM, "swarm_"));
  PetscCall(DMSetFromOptions(swarmDM));
  PetscCall(DMSetType(swarmDM, DMSWARM));
  PetscCall(PetscObjectSetName((PetscObject)swarmDM, "Swarm"));
  // Synchronize the swarm DM with the fluid DM.
  PetscCall(DMGetDimension(ctx->fluidDM, &dim));
  PetscCall(DMSetDimension(swarmDM, dim));
  PetscCall(DMSwarmSetCellDM(swarmDM, ctx->fluidDM));
  // Declare this to be a PIC swarm. This must occur after setting `dim`.
  PetscCall(DMSwarmSetType(swarmDM, DMSWARM_PIC));
  // Register non-default fields that each particle will have.
  PetscCall(DMSwarmInitializeFieldRegister(swarmDM));
  // --> (x, y, z) velocity components
  PetscCall(DMSwarmRegisterPetscDatatypeField(swarmDM, "velocity", ndim, PETSC_REAL));
  PetscCall(DMSwarmFinalizeFieldRegister(swarmDM));
  // Set the local number of points in each dimension.
  PetscCall(DMDAGetCorners(ctx->fluidDM, NULL, NULL, NULL, &ctx->grid.nx, &ctx->grid.ny, &ctx->grid.nz));
  // Set the number of charged particles equal to the default, if necessary.
  if (ctx->plasma.Np == -1) {
    ctx->plasma.Np = ctx->grid.Nx * ctx->grid.Ny * ctx->grid.Nz;
  }
  // Set the per-processor swarm size and buffer length for efficient resizing.
  np = (PetscInt)(ctx->plasma.Np / ctx->mpi.size);
  bufsize = (PetscInt)(0.25 * np);
  PetscCall(DMSwarmSetLocalSizes(swarmDM, np, bufsize));
  // View information about the swarm DM.
  {
    /* NOTE: This is work-around for the fact that there appears to be no way to
    request -*_dm_view for DMSwarm objects from the command line. See
    ${PETSC_DIR}/src/dm/impls/swarm/swarm.c::DMInitialize_Swarm */
    PetscBool requested, found;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-swarm_dm_view", &requested, &found));
    if (found && requested) {
      PetscCall(DMView(swarmDM, PETSC_VIEWER_STDOUT_WORLD));
    }
  }
  // Assign the swarm DM to the application context.
  ctx->swarmDM = swarmDM;

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}

