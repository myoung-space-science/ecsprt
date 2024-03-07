#include <stdarg.h>
#include <petsc.h>

PetscErrorCode printNone(const char *format, ...)
{
  PetscFunctionBeginUser;
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode printWorld(const char *format, ...)
{
  va_list args;
  char    message[4096]="";

  PetscFunctionBeginUser;

  va_start(args, format);
  vsprintf(message, format, args);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s", message));
  va_end(args);

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode printSelf(const char *format, ...)
{
  va_list args;
  char    message[4096]="";

  PetscFunctionBeginUser;

  va_start(args, format);
  vsprintf(message, format, args);
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "%s", message));
  va_end(args);

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode printRanks(const char *format, ...)
{
  va_list args;
  char    message[4096]="";

  PetscFunctionBeginUser;

  va_start(args, format);
  vsprintf(message, format, args);
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", message));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
  va_end(args);

  PetscFunctionReturn(PETSC_SUCCESS);
}

