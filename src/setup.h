#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <petsc.h>
#include "hybrid.h"

extern PetscErrorCode SetUpContext(CLI cli, Context *ctx);
extern PetscErrorCode DestroyContext(Context *ctx);
extern PetscErrorCode CreatePotentialDM(Context *ctx);
extern PetscErrorCode CreateGridDM(Context *ctx);
extern PetscErrorCode CreateIonsDM(Context *ctx);

#endif // INITIALIZE_H