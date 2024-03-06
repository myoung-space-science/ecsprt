#include <time.h>
#include <unistd.h>
#include <petsc.h>
#include "hybrid.h"

PetscErrorCode initialize(int argc, char **args, const char *help, Context *ctx)
{

  /* Initialize PETSc and MPI. */
  if (access("petsc.ini", F_OK) != -1) {
    PetscCall(PetscInitialize(&argc, &args, "petsc.ini", help));
  } else {
    PetscCall(PetscInitialize(&argc, &args, NULL, help));
  }
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &ctx->mpi.rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &ctx->mpi.size));

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode finalize(time_t startTime, time_t endTime)
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
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n----------------------------------------\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Start time: %s", asctime(localtime(&startTime))));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "End time:   %s", asctime(localtime(&endTime))));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Elapsed time: %4.1f %s\n", elapsedTime, timeUnit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------\n"));

  /* Finalize PETSc and MPI. */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n***************** END ******************\n"));
  PetscCall(PetscFinalize());

  PetscFunctionReturn(PETSC_SUCCESS);
}