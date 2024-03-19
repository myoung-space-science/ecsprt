/* 3-D Electrostatic Potential Solver */
static char help[] = "A tool for solving the 3D quasineutral electrostatic-potential equation.";

#include <time.h>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>
#include <slepceps.h>
#include "ecsprt.h"
#include "common.h"
#include "parameters.h"
#include "setup.h"
#include "potential.h"
#include "file-io.h"
#include "lhs.h"


typedef struct {
  char      inpath[PETSC_MAX_PATH_LEN]; // path to input file
  PetscReal fluxScale[NDIM];            // values by which to scale density into flux
  PetscBool viewLHS;                    // option to view LHS operator structure
} Application;


/* Process command-line arguments specific to the potential solver. */
PetscErrorCode ProcessSolverOptions(Context ctx, Application *app)
{
  char      strArg[PETSC_MAX_PATH_LEN]="";
  PetscReal realArg;
  PetscBool boolArg;
  PetscBool found;

  PetscFunctionBeginUser;
  ctx.log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(PetscOptionsGetString(NULL, NULL, "--input", strArg, sizeof(strArg), &found));
  if (found) {
    PetscCall(PetscStrcat(app->inpath, strArg));
  }
  PetscCall(PetscOptionsGetBool(NULL, NULL, "--view-lhs", &boolArg, &found));
  if (found) {
    app->viewLHS = boolArg;
  } else {
    app->viewLHS = PETSC_FALSE;
  }
  PetscCall(PetscOptionsGetReal(NULL, NULL, "--flux-scale", &realArg, &found));
  if (found) {
    app->fluxScale[0] = realArg;
    app->fluxScale[1] = realArg;
    app->fluxScale[2] = realArg;
  } else {
    app->fluxScale[0] = 1.0;
    app->fluxScale[1] = 1.0;
    app->fluxScale[2] = 1.0;
  }
  PetscCall(PetscOptionsGetReal(NULL, NULL, "--x-flux-scale", &realArg, &found));
  if (found) {
    app->fluxScale[0] = realArg;
  }
  PetscCall(PetscOptionsGetReal(NULL, NULL, "--y-flux-scale", &realArg, &found));
  if (found) {
    app->fluxScale[1] = realArg;
  }
  PetscCall(PetscOptionsGetReal(NULL, NULL, "--z-flux-scale", &realArg, &found));
  if (found) {
    app->fluxScale[2] = realArg;
  }

  ctx.log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode EchoSetup(Context ctx, Application app)
{
  PetscViewer viewer;

  PetscFunctionBeginUser;
  ctx.log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Echo common parameter values.
  PetscCall(EchoOptions(ctx));

  // Open the viewer in "append" mode.
  PetscCall(OpenASCIIAppend(PETSC_COMM_WORLD, ctx.optionsLog, &viewer, &ctx));

  // View simulation-specific parameter values.
  PetscCall(PetscViewerASCIIPrintf(viewer, "\n\n#Application-Specific Parameter Values\n"));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "-------------------------------------\n"));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "x flux scale = %f\n", app.fluxScale[0]));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "y flux scale = %f\n", app.fluxScale[1]));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "z flux scale = %f\n", app.fluxScale[2]));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "density input = %s\n", app.inpath));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "#End of Application-Specific Parameter Values\n"));

  // Free memory.
  PetscCall(PetscViewerDestroy(&viewer));

  ctx.log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


int main(int argc, char **args)
{
  time_t        startTime, endTime;
  CLI           cli;
  Application   app;
  Context       ctx;
  KSP           ksp;

  PetscFunctionBeginUser;

  /* Log start time. */
  time(&startTime);

  /* Initialize PETSc and MPI. */
  initialize(argc, args, help, &ctx);

  /* Initialize SLEPc. */
  PetscCall(SlepcInitialize(&argc, &args, (char *)0, help));

  /* Assign parameter values from user arguments or defaults. */
  ctx.log.status("Processing common options\n");
  PetscCall(ProcessOptions(&cli));

  /* Set up the common application context. */
  ctx.log.status("Setting up common parameters\n");
  PetscCall(SetUpContext(cli, &ctx));

  /* Process application-specific options. */
  ctx.log.status("Processing application-specific options\n");
  PetscCall(ProcessSolverOptions(ctx, &app));

  /* Echo the initial state. */
  ctx.log.status("Echoing parameter values to %s\n", ctx.optionsLog);
  PetscCall(EchoSetup(ctx, app));

  /* Set up the fluid grid. */
  ctx.log.status("Creating the fluid-grid DM\n");
  PetscCall(CreateGridDM(cli.ndim, &ctx));

  /* Read density and fluxes from disk. */
  ctx.log.status("Loading fluid quantities\n");
  PetscCall(LoadFluidQuantities(app.fluxScale, app.inpath, &ctx));

  /* Set up the discrete grid for the electrostatic potential. */
  ctx.log.status("Creating the electrostatic-potential DM\n");
  PetscCall(CreatePotentialDM(cli.ndim, &ctx));

  /* Set up the Krylov-solver context for the electrostatic potential. */
  ctx.log.status("Setting up the potential solver\n");
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetDM(ksp, ctx.potential.dm));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetComputeInitialGuess(ksp, ComputeInitialPhi, &ctx));
  PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, &ctx));
  PetscCall(KSPSetComputeOperators(ksp, ComputeLHS, &ctx));

  /* Compute the electrostatic potential. */
  ctx.log.status("Computing initial potential\n");
  PetscCall(ComputePotential(ksp, &ctx));

  /* Output arrays. */
  ctx.log.status("Writing fluid quantities to HDF5\n");
  PetscCall(OutputFluidHDF5("", &ctx));

  /* View the LHS matrix structure. */
  {
    PetscBool found, true;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "--view-lhs", &true, &found));
    if (found && true) {
      ctx.log.status("Writing LHS structure to binary\n");
      PetscCall(ViewLHS(ksp, &ctx));
    }
  }

  /* Compute LHS eigenvalues. */
  {
    PetscBool found, true;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "--lhs-eigenvalues", &true, &found));
    if (found && true) {
    ctx.log.status("Computing LHS eigenvalues\n");
      PetscCall(ComputeLHSEigenvalues(ksp, &ctx));
    }
  }

  /* Free memory. */
  ctx.log.status("Freeing objects\n");
  PetscCall(KSPDestroy(&ksp));
  PetscCall(DestroyContext(&ctx));

  /* Finalize SLEPc. */
  PetscCall(SlepcFinalize());

  /* Log end time. */
  time(&endTime);

  /* Complete final tasks. */
  finalize(startTime, endTime, &ctx);

  return 0;
}
