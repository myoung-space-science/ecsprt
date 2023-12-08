#include <unistd.h>
#include <petsc.h>
#include "hybrid.h"

void initialize(int argc, char **args, const char *help, Context *ctx)
{

  /* Initialize PETSc and MPI. */
  if (access("petsc.ini", F_OK) != -1) {
    PetscCall(PetscInitialize(&argc, &args, "petsc.ini", help));
  } else {
    PetscCall(PetscInitialize(&argc, &args, NULL, help));
  }
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &ctx->mpi.rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &ctx->mpi.size));

}

