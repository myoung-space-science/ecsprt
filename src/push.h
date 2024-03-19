#ifndef PUSH_H
#define PUSH_H

#include "ecsprt.h"

extern PetscErrorCode BorisMoverBB(PetscReal dt, Context *ctx);
extern PetscErrorCode BorisMoverBz(PetscReal dt, Context *ctx);

#endif // PUSH_H