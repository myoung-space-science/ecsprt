#include <time.h>
#include <unistd.h>
#include <petsc.h>
#include "ecsprt.h"

PetscErrorCode Initialize(int argc, char **args, const char *help, const char *name, Context *ctx)
{
  PetscInt  logLevel=1;
  PetscBool requested;
  PetscBool found;

  /* Initialize PETSc and MPI. */
  if ((argc > 1) && (access(args[1], F_OK) != -1)) {
    PetscCall(PetscInitialize(&argc, &args, args[1], help));
  } else if (access("petsc.ini", F_OK) != -1) {
    PetscCall(PetscInitialize(&argc, &args, "petsc.ini", help));
  } else {
    PetscCall(PetscInitialize(&argc, &args, NULL, help));
  }
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &ctx->mpi.rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &ctx->mpi.size));

  /* Echo version number and exit, if requested. */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "--version", &requested, &found));
  if (found) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s: %s\n", ACRONYM, PROJECT));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "[%s] Version %s\n\n", name, VERSION));
    PetscCall(PetscEnd());
  }

  /* Assign the logging functions. */
  ctx->log.world        = printNone;
  ctx->log.self         = printNone;
  ctx->log.ranks        = printNone;
  ctx->log.checkpoint   = printNone;
  ctx->log.status       = printNone;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "--log-level", &logLevel, &found));
  if (logLevel > 0) {
    ctx->log.status     = printWorld;
  }
  if (logLevel > 1) {
    ctx->log.world      = printWorld;
    ctx->log.self       = printSelf;
    ctx->log.ranks      = printRanks;
  }
  if (logLevel > 2) {
    ctx->log.checkpoint = printWorld;
  }

  ctx->log.status("[%s] Running on %d MPI process(es)\n\n", name, ctx->mpi.size);
  ctx->log.status("\n**************** START *****************\n\n");

  PetscFunctionReturn(PETSC_SUCCESS);
}


int Finalize(time_t startTime, time_t endTime, Context *ctx)
{
  float       elapsedTime;
  char        timeUnit[16]="";

  /* Log end time. */
  elapsedTime = (float)(endTime-startTime);
  if (elapsedTime >= 86400.0) {
    elapsedTime /= 86400.0;
    PetscCall(PetscStrcat(timeUnit, "day(s)"));
  } else if (elapsedTime >= 3600.0) {
    elapsedTime /= 3600.0;
    PetscCall(PetscStrcat(timeUnit, "hour(s)"));
  } else if (elapsedTime >= 60.0) {
    elapsedTime /= 60.0;
    PetscCall(PetscStrcat(timeUnit, "minute(s)"));
  } else {
    PetscCall(PetscStrcat(timeUnit, "second(s)"));
  }
  ctx->log.status("\n----------------------------------------\n");
  ctx->log.status("Start time: %s", asctime(localtime(&startTime)));
  ctx->log.status("End time:   %s", asctime(localtime(&endTime)));
  ctx->log.status("Elapsed time: %4.1f %s\n", elapsedTime, timeUnit);
  ctx->log.status("----------------------------------------\n");

  /* Finalize PETSc and MPI. */
  ctx->log.status("\n***************** END ******************\n");
  PetscCall(PetscFinalize());

  return 0;
}