#ifndef LOGGING_H
#define LOGGING_H

#include <petsc.h>

extern PetscErrorCode printNone(const char *message);
extern PetscErrorCode printWorld(const char *message);
extern PetscErrorCode printSelf(const char *message);
extern PetscErrorCode printRanks(const char *message);

typedef PetscErrorCode (*LogFunction)(const char *message);

#endif // LOGGING_H
