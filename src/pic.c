/* 3D Hybrid PIC */
static char help[] = "A 3D hybrid particle-in-cell (PIC) simulation.";

#include <time.h>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>
#include "ecsprt.h"
#include "common.h"
#include "parameters.h"
#include "setup.h"
#include "positions.h"
#include "velocities.h"
#include "potential.h"
#include "file-io.h"
#include "pic.h"

typedef struct {
  PetscInt    Nt;          // number of time steps
  PetscReal   dt;          // time-step width
  PetscInt    it;          // time-step counter
  PetscInt    Dt;          // output cadence
  PDistType PDistType; // type of initial position distribution
  VDistType VDistType; // type of initial velocity distribution
} Application;

/* Process command-line arguments specific to the PIC simulation. */
PetscErrorCode ProcessPICOptions(Context ctx, Application *app)
{
  PetscReal realArg;
  PetscEnum enumArg;
  PetscInt  intArg;
  PetscBool found;

  PetscFunctionBeginUser;
  ctx.log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(PetscOptionsGetEnum(NULL, NULL, "--position-dist", PDistTypes, &enumArg, &found));
  if (found) {
    app->PDistType = enumArg;
  } else {
    app->PDistType = PDIST_FLAT_SOBOL;
  }
  PetscCall(PetscOptionsGetEnum(NULL, NULL, "--velocity-dist", VDistTypes, &enumArg, &found));
  if (found) {
    app->VDistType = enumArg;
  } else {
    app->VDistType = VDIST_NORMAL;
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
    app->dt = 1.0 / ctx.ions.nu;
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
  PetscCall(PetscViewerASCIIPrintf(viewer,     "Nt = %d\n", app.Nt));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "dt = %f [s]\n", app.dt));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "initial positions = %s\n", PDistTypes[app.PDistType]));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "initial velocities = %s\n", VDistTypes[app.VDistType]));
  PetscCall(PetscViewerASCIIPrintf(viewer,     "#End of Application-Specific Parameter Values\n"));

  // Free memory.
  PetscCall(PetscViewerDestroy(&viewer));

  ctx.log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
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

  /* Log start time. */
  time(&startTime);

  /* Initialize PETSc and MPI. */
  initialize(argc, args, help, &ctx);

  /* Assign parameter values from user arguments or defaults. */
  ctx.log.status("Processing common options\n");
  PetscCall(ProcessOptions(&cli));

  /* Set up the common application context. */
  ctx.log.status("Setting up common parameters\n");
  PetscCall(SetUpContext(cli, &ctx));

  /* Process application-specific options. */
  ctx.log.status("Processing application-specific options\n");
  PetscCall(ProcessPICOptions(ctx, &app));

  /* Echo the initial state. */
  ctx.log.status("Echoing parameter values to %s\n", ctx.optionsLog);
  PetscCall(EchoSetup(ctx, app));

  /* Echo this stage. */
  ctx.log.status("\n=== Initial stage ===\n\n");

  /* Set up the fluid grid. */
  ctx.log.status("Creating the fluid-grid DM\n");
  PetscCall(CreateGridDM(cli.ndim, &ctx));

  /* Set up the ion swarm. */
  ctx.log.status("Creating the particle-swarm DM\n");
  PetscCall(CreateIonsDM(cli.ndim, &ctx));

  /* Set initial particle positions. */
  ctx.log.status("Initializing positions\n");
  PetscCall(InitializePositions(cli.ndim, app.PDistType, &ctx));

  /* Set initial particle velocities. */
  ctx.log.status("Initializing velocities\n");
  PetscCall(InitializeVelocities(cli.ndim, app.VDistType, &ctx));

  /* Compute initial density and flux.*/
  ctx.log.status("Collecting initial moments\n");
  PetscCall(ctx.ions.collect(&ctx));

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

  /* Compute initial electrostatic potential. */
  ctx.log.status("Computing initial potential\n");
  PetscCall(ComputePotential(ksp, &ctx));

  /* Create a format string for the time step. */
  itwidth = 1+PetscLog10Real(app.Nt);
  sprintf(itfmt, "%%0%dd", itwidth);

  /* Create a template for time-dependent filenames. */
  PetscCall(PetscStrcat(pathfmt, "-"));
  PetscCall(PetscStrcat(pathfmt, itfmt));

  /* Output initial conditions. */
  ctx.log.status("Writing initial fluid quantities to HDF5\n");
  PetscCall(OutputFluidHDF5("-initial", &ctx));
  ctx.log.status("Writing initial particle quantities to binary\n");
  PetscCall(OutputSwarmBinary("-initial", &ctx));

  /* Create a template for the time-step string. */
  PetscCall(PetscStrcat(stepfmt, "< Time step "));
  PetscCall(PetscStrcat(stepfmt, itfmt));
  PetscCall(PetscStrcat(stepfmt, " >\n"));

  /* Echo this stage. */
  ctx.log.status("\n\n=== Main time-step loop ===\n\n");
  /* Begin main time-step loop. */
  for (it=0; it<app.Nt; it++) {

    /* Create a string to display time step with the appropriate width. */
    sprintf(stepstr, stepfmt, it);
    ctx.log.status("%s", stepstr);

    /* Update velocities. */
    ctx.log.status("Updating velocities\n");
    PetscCall(UpdateVelocities(cli.ndim, app.dt, &ctx));

    /* Update positions. */
    ctx.log.status("Updating positions\n");
    PetscCall(UpdatePositions(cli.ndim, app.dt, &ctx));

    /* Compute density and flux from ion positions. */
    ctx.log.status("Collecting moments\n");
    PetscCall(ctx.ions.collect(&ctx));

    /* Compute potential from density. */
    ctx.log.status("Computing potential\n");
    PetscCall(ComputePotential(ksp, &ctx));

    /* Output current time step. */
    if ((app.Dt > 0) && (it % app.Dt == 0)) {
      sprintf(pathstr, pathfmt, it);
      ctx.log.status("Writing fluid quantities to HDF5\n");
      PetscCall(OutputFluidHDF5(pathstr, &ctx));
      ctx.log.status("Writing particle quantities to binary\n");
      PetscCall(OutputSwarmBinary(pathstr, &ctx));
    }

    /* Print a newline to separate this time step from the next. */
    ctx.log.status("\n");
  }

  /* Echo this stage. */
  ctx.log.status("\n=== Final stage ===\n\n");

  /* Output final conditions. */
  ctx.log.status("Writing final fluid quantities to HDF5\n");
  PetscCall(OutputFluidHDF5("-final", &ctx));
  ctx.log.status("Writing final particle quantities to binary\n");
  PetscCall(OutputSwarmBinary("-final", &ctx));

  /* Free memory. */
  ctx.log.status("Freeing objects\n");
  PetscCall(KSPDestroy(&ksp));
  PetscCall(DestroyContext(&ctx));
  PetscCall(DMDestroy(&ctx.swarmDM));

  /* Log end time. */
  time(&endTime);

  /* Complete final tasks. */
  finalize(startTime, endTime, &ctx);

  return 0;
}