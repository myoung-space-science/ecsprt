#include <petsc.h>
#include "hybrid.h"


/* Set parameter values common to simulation and solver applications. */
PetscErrorCode SetUpContext(CLI cli, Context *ctx)
{

  PetscReal tmp;
  PetscInt  ib;

  PetscFunctionBeginUser;
  ECHO_FUNCTION_ENTER;

  // Declare the name of the options log.
  PetscCall(PetscStrcpy(ctx->optionsLog, "options.log"));

  // Set fundamental parameter values.
  ctx->electrons.q = -Q;
  ctx->electrons.m = ME;

  // Explicitly set the neutral charge to zero.
  ctx->neutrals.q = 0.0;

  // Set the LHS function based on LHS type.
  switch (cli.lhsType) {
  case LHS_IDENTITY:
    ctx->potential.lhs = ComputeIdentityLHS;
    ctx->potential.stencilSize = 1;
    ctx->potential.stencilType = DMDA_STENCIL_STAR;
    break;
  case LHS_LAPLACIAN:
    ctx->potential.lhs = ComputeLaplacianLHS;
    ctx->potential.stencilSize = 7;
    ctx->potential.stencilType = DMDA_STENCIL_STAR;
    break;
  case LHS_FULL:
    ctx->potential.lhs = ComputeFullLHS;
    ctx->potential.stencilSize = 11;
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
    ctx->potential.rhs = ComputeSinusoidalRHS;
    break;
  case RHS_FULL:
    ctx->potential.rhs = ComputeFullRHS;
    break;
  default:
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown RHS type: \"%s\"\n", RHSTypes[cli.rhsType]);
  }

  // Copy plasma and species parameter values.
  ctx->plasma.n0 = cli.n0;
  ctx->plasma.B0 = cli.B0;
  ctx->plasma.E0 = cli.E0;
  ctx->plasma.Np = cli.Np;
  ctx->ions.q = cli.qi;
  ctx->ions.m = cli.mi;
  ctx->ions.v0x = cli.vi0x;
  ctx->ions.v0y = cli.vi0y;
  ctx->ions.v0z = cli.vi0z;
  ctx->ions.v0 = cli.vi0;
  ctx->ions.vTx = cli.viTx;
  ctx->ions.vTy = cli.viTy;
  ctx->ions.vTz = cli.viTz;
  ctx->ions.vT = cli.viT;
  ctx->ions.T = cli.Ti;
  ctx->ions.nu = cli.nui;
  ctx->electrons.v0x = cli.ve0x;
  ctx->electrons.v0y = cli.ve0y;
  ctx->electrons.v0z = cli.ve0z;
  ctx->electrons.v0 = cli.ve0;
  ctx->electrons.vTx = cli.veTx;
  ctx->electrons.vTy = cli.veTy;
  ctx->electrons.vTz = cli.veTz;
  ctx->electrons.vT = cli.veT;
  ctx->electrons.T = cli.Te;
  ctx->electrons.nu = cli.nue;
  ctx->neutrals.m = cli.mn;
  ctx->neutrals.v0x = cli.vn0x;
  ctx->neutrals.v0y = cli.vn0y;
  ctx->neutrals.v0z = cli.vn0z;
  ctx->neutrals.v0 = cli.vn0;
  ctx->neutrals.vTx = cli.vnTx;
  ctx->neutrals.vTy = cli.vnTy;
  ctx->neutrals.vTz = cli.vnTz;
  ctx->neutrals.vT = cli.vnT;
  ctx->neutrals.T = cli.Tn;

  // Copy grid parameters.
  ctx->grid.Nx = cli.Nx;
  ctx->grid.Ny = cli.Ny;
  ctx->grid.Nz = cli.Nz;
  ctx->grid.dx = cli.dx;
  ctx->grid.dy = cli.dy;
  ctx->grid.dz = cli.dz;
  ctx->grid.x0 = cli.x0;
  ctx->grid.y0 = cli.y0;
  ctx->grid.z0 = cli.z0;
  ctx->grid.x1 = cli.x1;
  ctx->grid.y1 = cli.y1;
  ctx->grid.z1 = cli.z1;

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
  if ((cli.xBT[0] == BT_PERIODIC) && (cli.xBT[1] == BT_PERIODIC)) {
    ctx->ions.xBC = BC_PERIODIC;
    ctx->grid.xBC = DM_BOUNDARY_PERIODIC;
  } else if ((cli.xBT[0] == BT_INJECTION) && (cli.xBT[1] == BT_REFLECTION)) {
    ctx->ions.xBC = BC_INJECT_REFLECT;
    ctx->grid.xBC = DM_BOUNDARY_GHOSTED;
  } else if ((cli.xBT[0] == BT_INJECTION) && (cli.xBT[1] == BT_ADVECTION)) {
    ctx->ions.xBC = BC_INJECT_ADVECT;
    ctx->grid.xBC = DM_BOUNDARY_GHOSTED;
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP, "Inconsistent x-dimension boundary conditions: {%s, %s}\n", BoundaryTypes[cli.xBT[0]], BoundaryTypes[cli.xBT[1]]);
  }
  // y dimension
  if ((cli.yBT[0] == BT_PERIODIC) && (cli.yBT[1] == BT_PERIODIC)) {
    ctx->ions.yBC = BC_PERIODIC;
    ctx->grid.yBC = DM_BOUNDARY_PERIODIC;
  } else if ((cli.yBT[0] == BT_INJECTION) && (cli.yBT[1] == BT_REFLECTION)) {
    ctx->ions.yBC = BC_INJECT_REFLECT;
    ctx->grid.yBC = DM_BOUNDARY_GHOSTED;
  } else if ((cli.yBT[0] == BT_INJECTION) && (cli.yBT[1] == BT_ADVECTION)) {
    ctx->ions.yBC = BC_INJECT_ADVECT;
    ctx->grid.yBC = DM_BOUNDARY_GHOSTED;
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP, "Inconsistent y-dimension boundary conditions: {%s, %s}\n", BoundaryTypes[cli.yBT[0]], BoundaryTypes[cli.yBT[1]]);
  }
  // z dimension
  if ((cli.zBT[0] == BT_PERIODIC) && (cli.zBT[1] == BT_PERIODIC)) {
    ctx->ions.zBC = BC_PERIODIC;
    ctx->grid.zBC = DM_BOUNDARY_PERIODIC;
  } else if ((cli.zBT[0] == BT_INJECTION) && (cli.zBT[1] == BT_REFLECTION)) {
    ctx->ions.zBC = BC_INJECT_REFLECT;
    ctx->grid.zBC = DM_BOUNDARY_GHOSTED;
  } else if ((cli.zBT[0] == BT_INJECTION) && (cli.zBT[1] == BT_ADVECTION)) {
    ctx->ions.zBC = BC_INJECT_ADVECT;
    ctx->grid.zBC = DM_BOUNDARY_GHOSTED;
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP, "Inconsistent z-dimension boundary conditions: {%s, %s}\n", BoundaryTypes[cli.zBT[0]], BoundaryTypes[cli.zBT[1]]);
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
    PRINT_WORLD("Warning: Setting neutral temperature equal to ion temperature (%.1f K)\n", ctx->ions.T);
  }

  // Set neutral thermal velocity from temperature.
  ctx->neutrals.vT = PetscSqrtReal(2.0 * KB * ctx->neutrals.T / ctx->neutrals.m);

  // [DEV] Hard-code the electron thermal coefficient to 1.0 (isothermal).
  ctx->gammaT = 1.0;

  // TODO: Should we set default collision frequencies based on an analytic
  // formulation (e.g., from Schunk & Nagy)?

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Destroy the application context and associated objects. */
PetscErrorCode DestroyContext(Context *ctx)
{

  PetscFunctionBeginUser;
  ECHO_FUNCTION_ENTER;

  PetscCall(VecDestroy(&ctx->moments));
  PetscCall(DMDestroy(&ctx->swarmDM));

  ECHO_FUNCTION_EXIT;
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
PetscErrorCode CreateGridDM(Context *ctx)
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
  ECHO_FUNCTION_ENTER;

  // Create the DM.
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, xBC, yBC, zBC, stencilType, Nx, Ny, Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, NULL, &dm));
  // Perform basic setup.
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)dm, "grid_"));
  PetscCall(DMDASetElementType(dm, DMDA_ELEMENT_Q1));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(PetscObjectSetName((PetscObject)dm, "GridDM"));
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
  // Set the local number of points in each dimension.
  PetscCall(DMDAGetCorners(dm, NULL, NULL, NULL, &ctx->grid.nx, &ctx->grid.ny, &ctx->grid.nz));
  // Set the number of charged particles equal to the default, if necessary.
  if (ctx->plasma.Np == -1) {
    ctx->plasma.Np = ctx->grid.Nx * ctx->grid.Ny * ctx->grid.Nz;
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

  // Set uniform coordinates on the DM.
  PetscCall(DMDASetUniformCoordinates(dm, ctx->grid.x0, ctx->grid.x1, ctx->grid.y0, ctx->grid.y1, ctx->grid.z0, ctx->grid.z1));
  // Declare grid-quantity names.
  PetscCall(DMDASetFieldName(dm, 0, "density"));
  PetscCall(DMDASetFieldName(dm, 1, "x flux"));
  PetscCall(DMDASetFieldName(dm, 2, "y flux"));
  PetscCall(DMDASetFieldName(dm, 3, "z flux"));
  // View information about the DM.
  PetscCall(DMView(dm, PETSC_VIEWER_STDOUT_WORLD));
  // Create a persistent vector for outputing fluid quantities.
  PetscCall(DMCreateGlobalVector(dm, &ctx->moments));
  PetscCall(VecZeroEntries(ctx->moments));
  // Assign the grid DM to the application context.
  ctx->fluidDM = dm;

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Create the data manager for the electrostatic potential.

The logical grid for the electrostatic potential will have the same shape,
stencil type, and boundary conditions as the logical grid for fluid quantities.
Unlike the latter, it will use a stencil width of 2 in order to accommodate
second-order finite differencing of the electric field.

This function should follow `CreateGridDM`, which performs various set-up and
sychronization tasks on grid parameters.
*/
PetscErrorCode CreatePotentialDM(Context *ctx)
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
  ECHO_FUNCTION_ENTER;

  // Create the DM object.
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, xBC, yBC, zBC, stencilType, Nx, Ny, Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, NULL, &ctx->potential.dm));
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
  // Echo information about the DM.
  PetscCall(DMView(ctx->potential.dm, PETSC_VIEWER_STDOUT_WORLD));

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Create the data manager for ion distributions.

The logical grid for ion distributions unconditionally uses ghosted boundaries
(i.e., regardless of the user-requested boundary types). Downstream functions
are responsible for applying the appropriate boundary conditions.
*/
PetscErrorCode CreateIonsDM(Context *ctx)
{
  PetscInt        Nx=ctx->grid.Nx;
  PetscInt        Ny=ctx->grid.Ny;
  PetscInt        Nz=ctx->grid.Nz;
  DMBoundaryType  xBC=ctx->grid.xBC;
  DMBoundaryType  yBC=ctx->grid.yBC;
  DMBoundaryType  zBC=ctx->grid.zBC;
  DMDAStencilType stencilType=DMDA_STENCIL_BOX;
  PetscInt        dof=4;
  PetscInt        width=1;
  PetscInt        dim;
  DM              ionsDM, cellDM;
  PetscInt        bufsize=0;
  PetscInt        np;
  PetscReal       x0=ctx->grid.x0;
  PetscReal       y0=ctx->grid.y0;
  PetscReal       z0=ctx->grid.z0;
  PetscReal       x1=ctx->grid.x1+(2*ctx->grid.dx);
  PetscReal       y1=ctx->grid.y1+(2*ctx->grid.dy);
  PetscReal       z1=ctx->grid.z1+(2*ctx->grid.dz);

  PetscFunctionBeginUser;
  ECHO_FUNCTION_ENTER;

  // Create the cell DM.
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, xBC, yBC, zBC, stencilType, Nx, Ny, Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof, width, NULL, NULL, NULL, &cellDM));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)cellDM, "cell_"));
  PetscCall(DMDASetElementType(cellDM, DMDA_ELEMENT_Q1));
  PetscCall(DMSetFromOptions(cellDM));
  PetscCall(DMSetUp(cellDM));
  PetscCall(PetscObjectSetName((PetscObject)cellDM, "CellDM"));
  // Set uniform coordinates on the DM.
  PetscCall(DMDASetUniformCoordinates(cellDM, x0, x1, y0, y1, z0, z1));
  // Create the ions DM.
  PetscCall(DMCreate(PETSC_COMM_WORLD, &ionsDM));
  // Perform basic setup.
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)ionsDM, "ions_"));
  PetscCall(DMSetType(ionsDM, DMSWARM));
  PetscCall(PetscObjectSetName((PetscObject)ionsDM, "Ions"));
  // Synchronize the ions DM with the vlasov DM.
  PetscCall(DMGetDimension(cellDM, &dim));
  PetscCall(DMSetDimension(ionsDM, dim));
  PetscCall(DMSwarmSetCellDM(ionsDM, cellDM));
  // Declare this to be a PIC swarm. This must occur after setting `dim`.
  PetscCall(DMSwarmSetType(ionsDM, DMSWARM_PIC));
  // Register non-default fields that each particle will have.
  PetscCall(DMSwarmInitializeFieldRegister(ionsDM));
  // --> (x, y, z) velocity components
  PetscCall(DMSwarmRegisterPetscDatatypeField(ionsDM, "velocity", NDIM, PETSC_REAL));
  PetscCall(DMSwarmFinalizeFieldRegister(ionsDM));
  // Set the per-processor swarm size and buffer length for efficient resizing.
  np = (PetscInt)(ctx->plasma.Np / ctx->mpi.size);
  bufsize = (PetscInt)(0.25 * np);
  PetscCall(DMSwarmSetLocalSizes(ionsDM, np, bufsize));
  // View information about the cell DM.
  PetscCall(DMView(cellDM, PETSC_VIEWER_STDOUT_WORLD));
  // View information about the ions DM.
  PetscCall(DMView(ionsDM, PETSC_VIEWER_STDOUT_WORLD));
  // Assign the ions DM to the application context.
  ctx->swarmDM = ionsDM;

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}

