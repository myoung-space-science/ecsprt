#ifndef LOGGING_H
#define LOGGING_H

#include <petsc.h>

extern PetscErrorCode printNone(const char *format, ...);
extern PetscErrorCode printWorld(const char *format, ...);
extern PetscErrorCode printSelf(const char *format, ...);
extern PetscErrorCode printRanks(const char *format, ...);

typedef PetscErrorCode (*LogFunction)(const char *format, ...);

#endif // LOGGING_H
