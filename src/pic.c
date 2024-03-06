/* 3D Hybrid PIC */
static char help[] = "A 3D hybrid particle-in-cell (PIC) simulation.";

#include <time.h>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>
#include "hybrid.h"
#include "common.h"
#include "parameters.h"
#include "setup.h"
#include "distributions.h"
#include "particles.h"
#include "potential.h"
#include "file-io.h"
#include "pic.h"

typedef struct {
  PetscInt    Nt;          // number of time steps
  PetscReal   dt;          // time-step width
  PetscInt    it;          // time-step counter
  PetscInt    Dt;          // output cadence
  DensityType densityType; // type of initial density
} Application;

/* Process command-line arguments specific to the PIC simulation. */
PetscErrorCode ProcessPICOptions(Context ctx, Application *app)
{
  PetscReal realArg;
  PetscEnum enumArg;
  PetscInt  intArg;
  PetscBool found;

  PetscFunctionBeginUser;
  ECHO_FUNCTION_ENTER;

  PetscCall(PetscOptionsGetEnum(NULL, NULL, "--density-type", DensityTypes, &enumArg, &found));
  if (found) {
    app->densityType = enumArg;
  } else {
    app->densityType = DENSITY_FLAT_SOBOL;
  }
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-Nt", &intArg, &found));
  if (found) {
    app->Nt = intArg;
  } else {
    app->Nt = 1;
  }
  PetscCall(PetscOptionsGetInt(NULL, NULL, "--every", &intArg, &found));
  if (found) {
    app->Dt = intArg;
  } else {
    app->Dt = 1;
  }
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &realArg, &found));
  if (found) {
    app->dt = realArg;
  } else {
    PRINT_WORLD("Warning: Setting dt = 1 / nui\n");
    app->dt = 1.0 / ctx.ions.nu;
  }

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode EchoSetup(Context ctx, Application app)
{
  PetscViewer viewer;

  PetscFunctionBeginUser;
  ECHO_FUNCTION_ENTER;

  // Echo common parameter values.
  PetscCall(EchoOptions(ctx));

  // Open the viewer in "append" mode.
  PetscCall(OpenASCIIAppend(PETSC_COMM_WORLD, ctx.optionsLog, &viewer));

  // View simulation-specific parameter values.
  PetscCall(PetscViewerASCIIPrintf(viewer, "\n\n#Application-Specific Parameter Values\n"));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "Nt = %d\n", app.Nt));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "dt = %f [s]\n", app.dt));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "density type = %s\n", DensityTypes[app.densityType]));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "#End of Application-Specific Parameter Values\n"));

  // Free memory.
  PetscCall(PetscViewerDestroy(&viewer));

  ECHO_FUNCTION_EXIT;
  PetscFunctionReturn(PETSC_SUCCESS);
}


int main(int argc, char **args)
{
  time_t      startTime, endTime;
  CLI         cli;
  Application app;
  Context     ctx;
  KSP         ksp;
  PetscInt    it;
  char        itfmt[5];
  char        pathfmt[PETSC_MAX_PATH_LEN]="", pathstr[PETSC_MAX_PATH_LEN];
  char        stepfmt[256]="", stepstr[256];
  PetscInt    itwidth;

  PetscFunctionBeginUser;

  /* Initialize PETSc and MPI. */
  initialize(argc, args, help, &ctx);

  /* Assign parameter values from user arguments or defaults. */
  PetscCall(ProcessOptions(&cli));

  /* Log start time. */
  time(&startTime);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n**************** START *****************\n\n"));

  /* Echo this stage. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n*** Initialization stage ***\n\n"));

  /* Set up the common application context. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Processing common options\n"));
  PetscCall(SetUpContext(cli, &ctx));

  /* Process application-specific options. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Processing application-specific options\n"));
  PetscCall(ProcessPICOptions(ctx, &app));

  /* Set up the fluid grid. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Creating the fluid-grid DM\n"));
  PetscCall(CreateGridDM(&ctx));

  /* Set up the ion swarm. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Creating the particle-swarm DM\n"));
  PetscCall(CreateIonsDM(&ctx));

  /* Set initial particle positions. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Initializing positions\n"));
  PetscCall(InitializePositions(app.densityType, &ctx));

  /* Set initial particle velocities. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Initializing velocities\n"));
  PetscCall(InitializeVelocities(&ctx));

  /* Echo the initial state. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Echoing parameter values to %s\n", ctx.optionsLog));
  PetscCall(EchoSetup(ctx, app));

  /* Compute initial density and flux.*/
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Collecting initial moments\n"));
  PetscCall(CollectFluidMoments(&ctx));

  /* Set up the discrete grid for the electrostatic potential. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Creating the electrostatic-potential DM\n"));
  PetscCall(CreatePotentialDM(&ctx));

  /* Set up the Krylov-solver context for the electrostatic potential. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Setting up the potential solver\n"));
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetDM(ksp, ctx.potential.dm));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetComputeInitialGuess(ksp, ComputeInitialPhi, &ctx));
  PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, &ctx));
  PetscCall(KSPSetComputeOperators(ksp, ComputeLHS, &ctx));

  /* Compute initial electrostatic potential. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Computing initial potential\n"));
  PetscCall(ComputePotential(ksp, &ctx));

  /* Create a format string for the time step. */
  itwidth = 1+PetscLog10Real(app.Nt);
  sprintf(itfmt, "%%0%dd", itwidth);

  /* Create a template for time-dependent filenames. */
  PetscCall(PetscStrcat(pathfmt, "-"));
  PetscCall(PetscStrcat(pathfmt, itfmt));

  /* Output initial conditions. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing initial fluid quantities to HDF5\n"));
  PetscCall(OutputFluidHDF5("-initial", &ctx));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing initial particle quantities to binary\n"));
  PetscCall(OutputSwarmBinary("-initial", &ctx));

  /* Create a template for the time-step string. */
  PetscCall(PetscStrcat(stepfmt, "< Time step "));
  PetscCall(PetscStrcat(stepfmt, itfmt));
  PetscCall(PetscStrcat(stepfmt, " >\n"));

  /* Echo this stage. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n*** Main time-step loop ***\n\n"));
  /* Begin main time-step loop. */
  for (it=0; it<app.Nt; it++) {

    /* Create a string to display time step with the appropriate width. */
    sprintf(stepstr, stepfmt, it);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s", stepstr));

    /* Update velocities. */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Updating velocities\n"));
    PetscCall(UpdateVelocities(app.dt, &ctx));

    /* Update positions. */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Updating positions\n"));
    PetscCall(UpdatePositions(app.dt, &ctx));

    /* Compute density and flux from ion positions. */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Collecting moments\n"));
    PetscCall(CollectFluidMoments(&ctx));

    /* Compute potential from density. */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Computing potential\n"));
    PetscCall(ComputePotential(ksp, &ctx));

    /* Output current time step. */
    if ((app.Dt > 0) && (it % app.Dt == 0)) {
      sprintf(pathstr, pathfmt, it);
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing fluid quantities to HDF5\n"));
      PetscCall(OutputFluidHDF5(pathstr, &ctx));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing particle quantities to binary\n"));
      PetscCall(OutputSwarmBinary(pathstr, &ctx));
    }

    /* Print a newline to separate this time step from the next. */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

  /* Echo this stage. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n*** Finalization stage ***\n\n"));

  /* Output final conditions. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing final fluid quantities to HDF5\n"));
  PetscCall(OutputFluidHDF5("-final", &ctx));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing final particle quantities to binary\n"));
  PetscCall(OutputSwarmBinary("-final", &ctx));

  /* Free memory. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Freeing objects\n"));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(DestroyContext(&ctx));

  /* Log end time. */
  time(&endTime);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n----------------------------------------\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Start time: %s", asctime(localtime(&startTime))));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "End time:   %s", asctime(localtime(&endTime))));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Elapsed time: %f s\n", (float)(endTime-startTime)));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------\n"));

  /* Finalize PETSc and MPI. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n***************** END ******************\n"));
  PetscCall(PetscFinalize());

  return 0;
}