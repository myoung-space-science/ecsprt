#ifndef FILE_IO_H
#define FILE_IO_H

#include <petsc.h>
#include "hybrid.h"

extern PetscErrorCode LoadFluidQuantities(PetscReal fluxScale[NDIM], char inpath[PETSC_MAX_PATH_LEN], Context *ctx);
extern PetscErrorCode OpenASCIIAppend(MPI_Comm comm, const char filename[], PetscViewer *viewer);
extern PetscErrorCode OutputFluidHDF5(const char *insert, Context *ctx);
extern PetscErrorCode OutputSwarmBinary(const char *insert, Context *ctx);
extern PetscErrorCode ViewLHS(KSP ksp);

#endif // FILE_IO_H
