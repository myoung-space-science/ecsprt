#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <petsc.h>
#include "ecsprt.h"

extern PetscErrorCode SetUpContext(CLI cli, Context *ctx);
extern PetscErrorCode DestroyContext(Context *ctx);
extern PetscErrorCode CreateGridDM(PetscInt ndim, Context *ctx);
extern PetscErrorCode CreatePotentialDM(PetscInt ndim, Context *ctx);
extern PetscErrorCode CreateIonsDM(PetscInt ndim, Context *ctx);

#endif // INITIALIZE_H