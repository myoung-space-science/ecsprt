#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <petsc.h>
#include "ecsprt.h"

extern PetscErrorCode ProcessOptions(CLI *cli);
extern PetscErrorCode EchoOptions(Context ctx);

#endif // PARAMETERS_H
