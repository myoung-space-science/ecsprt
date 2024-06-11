#ifndef COMMON_H
#define COMMON_H

#include <petsc.h>
#include "ecsprt.h"

extern PetscErrorCode Initialize(int argc, char **args, const char *help, const char *name, Context *ctx);
extern int Finalize(Context *ctx);

#endif // COMMON_H
