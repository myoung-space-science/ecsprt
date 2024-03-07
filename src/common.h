#ifndef COMMON_H
#define COMMON_H

#include <petsc.h>
#include "hybrid.h"

extern PetscErrorCode initialize(int argc, char **args, const char *help, Context *ctx);
extern PetscErrorCode finalize(time_t startTime, time_t endTime, Context *ctx);

#endif // COMMON_H
