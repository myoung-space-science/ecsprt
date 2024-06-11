#ifndef FILE_IO_H
#define FILE_IO_H

#include <petsc.h>
#include "ecsprt.h"

extern PetscErrorCode LoadFluidQuantities(PetscReal fluxScale[3], char inpath[PETSC_MAX_PATH_LEN], Context *ctx);
extern PetscErrorCode OpenASCIIAppend(MPI_Comm comm, const char filename[], PetscViewer *viewer, Context *ctx);
extern PetscErrorCode OutputFluidHDF5(const char *insert, Context *ctx);
extern PetscErrorCode OutputSwarmBinary(const char *insert, Context *ctx);
extern PetscErrorCode OutputNoOp(const char *insert, Context *ctx);
extern PetscErrorCode ViewLHS(KSP ksp, Context *ctx);

typedef PetscErrorCode (*OutputFunc)(const char *s, Context *ctx);

#endif // FILE_IO_H
