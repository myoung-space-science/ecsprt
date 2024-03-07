/* 3-D Electrostatic Potential Solver */
static char help[] = "A tool for solving the 3D quasineutral electrostatic-potential equation.";

#include <time.h>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>
#include <slepceps.h>
#include "hybrid.h"
#include "common.h"
#include "parameters.h"
#include "setup.h"
#include "particles.h"
#include "potential.h"
#include "file-io.h"


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

  /* Initialize PETSc and MPI. */
  initialize(argc, args, help, &ctx);

  /* Assign parameter values from user arguments or defaults. */
  PetscCall(ProcessOptions(&cli));

  /* Log start time. */
  time(&startTime);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n**************** START *****************\n\n"));

  /* Initialize SLEPc. */
  PetscCall(SlepcInitialize(&argc, &args, (char *)0, help));

  /* Set up the common application context. */
  PetscCall(SetUpContext(cli, &ctx));

  /* Process application-specific options. */
  PetscCall(ProcessSolverOptions(ctx, &app));

  /* Set up the ions. */
  PetscCall(CreateIonsDM(&ctx));

  /* Echo the initial state. */
  PetscCall(EchoSetup(ctx, app));

  /* Read density and fluxes from disk. */
  PetscCall(LoadFluidQuantities(app.fluxScale, app.inpath, &ctx));

  /* Set up the discrete grid for the electrostatic potential. */
  PetscCall(CreatePotentialDM(&ctx));

  /* Set up the Krylov-solver context for the electrostatic potential. */
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetDM(ksp, ctx.potential.dm));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetComputeInitialGuess(ksp, ComputeInitialPhi, &ctx));
  PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, &ctx));
  PetscCall(KSPSetComputeOperators(ksp, ComputeLHS, &ctx));

  /* Compute the electrostatic potential. */
  PetscCall(ComputePotential(ksp, &ctx));

  /* Output arrays. */
  PetscCall(OutputFluidHDF5("", &ctx));

  /* View the LHS matrix structure. */
  {
    PetscBool found, true;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "--view-lhs", &true, &found));
    if (found && true) {
      PetscCall(ViewLHS(ksp, &ctx));
    }
  }

  /* Compute LHS eigenvalues. */
  {
    PetscBool found, true;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "--lhs-eigenvalues", &true, &found));
    if (found && true) {
      PetscCall(ComputeLHSEigenvalues(ksp, &ctx));
    }
  }

  /* Free memory. */
  PetscCall(KSPDestroy(&ksp));
  PetscCall(DestroyContext(&ctx));

  /* Log end time. */
  time(&endTime);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n----------------------------------------\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Start time: %s", asctime(localtime(&startTime))));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "End time:   %s", asctime(localtime(&endTime))));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Elapsed time: %f s\n", (float)(endTime-startTime)));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------\n"));

  /* Finalize SLEPc. */
  PetscCall(SlepcFinalize());

  /* Finalize PETSc and MPI. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n***************** END ******************\n"));
  PetscCall(PetscFinalize());

  return 0;
}
