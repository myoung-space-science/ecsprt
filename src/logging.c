#include <petsc.h>

PetscErrorCode printNone(const char *message)
{
  PetscFunctionBeginUser;
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode printWorld(const char *message)
{
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, message));
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode printSelf(const char *message)
{
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(PETSC_COMM_SELF, message));
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode printRanks(const char *message)
{
  PetscFunctionBeginUser;
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, message));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
  PetscFunctionReturn(PETSC_SUCCESS);
}

