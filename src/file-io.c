#include <petsc.h>
#include <petscviewerhdf5.h>
#include "hybrid.h"


PetscErrorCode LoadFluidQuantities(PetscReal fluxScale[NDIM], char inpath[PETSC_MAX_PATH_LEN], Context *ctx)
{
  PetscBool     nullPath;
  char          fullpath[PETSC_MAX_PATH_LEN];
  PetscViewer   viewer;
  DM           *dms, dm;
  DM            fluidDM=ctx->fluidDM;
  PetscInt      Nf;
  char        **keys;
  PetscInt      field;
  Vec           density, tmpflux, moments=ctx->moments;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Check for user-provided density file.
  PetscCall(PetscStrcmp(inpath, "", &nullPath));
  if (nullPath) {
    ctx->log.world("Warning: Got empty path to density file; using built-in sinusoidal form of RHS.\n");
    ctx->potential.rhs = ComputeSinusoidalRHS;
    PetscFunctionReturn(PETSC_SUCCESS);
  } else {
    PetscCall(PetscGetFullPath(inpath, fullpath, PETSC_MAX_PATH_LEN));
  }

  // Create the HDF5 viewer.
  PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, fullpath, FILE_MODE_READ, &viewer));

  // Zero the target vlasov vector.
  PetscCall(VecZeroEntries(moments));

  // Load density from the HDF5 file.
  PetscCall(DMCreateFieldDecomposition(fluidDM, &Nf, &keys, NULL, &dms));
  ctx->log.world("Attempting to load density from %s\n", fullpath);
  field = 0;
  dm = dms[field];
  PetscCall(DMGetGlobalVector(dm, &density));
  PetscCall(PetscObjectSetName((PetscObject)density, "density-kji"));
  PetscCall(VecLoad(density, viewer));
  PetscCall(VecStrideScatter(density, field, moments, INSERT_VALUES));
  ctx->log.world("Loaded density\n");

  // Convert density into fluxes.
  for (field=1; field<Nf; field++) {
    dm = dms[field];
    PetscCall(DMGetGlobalVector(dm, &tmpflux));
    PetscCall(VecZeroEntries(tmpflux));
    PetscCall(VecAXPY(tmpflux, fluxScale[field-1], density));
    PetscCall(VecStrideScatter(tmpflux, field, moments, INSERT_VALUES));
    ctx->log.world("Created %s from density\n", keys[field]);
    PetscCall(DMRestoreGlobalVector(dm, &tmpflux));
  }
  PetscCall(DMRestoreGlobalVector(dms[0], &density));

  // Release memory.
  for (field=0; field<Nf; field++) {
    PetscFree(keys[field]);
    PetscCall(DMDestroy(&dms[field]));
  }
  PetscFree(keys);
  PetscFree(dms);

  // Destroy the HDF5 viewer.
  PetscCall(PetscViewerDestroy(&viewer));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Open a PETSc ASCII viewer for appending text.

This is a convenience wrapper around standard PETSc routines. It is meant to
replace `PetscViewerASCIIOpen`, which opens the file for writting and clobbers
any existing contents.

As with `PetscViewerASCIIOpen`, you should call `PetscViewerDestroy` on the
viewer object when you no longer need it.
*/
PetscErrorCode OpenASCIIAppend(MPI_Comm comm, const char filename[], PetscViewer *viewer, Context *ctx)
{

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(PetscViewerCreate(comm, viewer));
  PetscCall(PetscViewerSetType(*viewer, PETSCVIEWERASCII));
  PetscCall(PetscViewerFileSetMode(*viewer, FILE_MODE_APPEND));
  PetscCall(PetscViewerFileSetName(*viewer, filename));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode OutputSwarmBinary(const char *insert, Context *ctx)
{
  char          posfn[PETSC_MAX_PATH_LEN]="position";
  char          velfn[PETSC_MAX_PATH_LEN]="velocity";
  PetscViewer   viewer;
  DM            swarmDM=ctx->swarmDM;
  Vec           target;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Do we need to call
  // https://petsc.org/release/manualpages/Viewer/PetscViewerBinarySetUseMPIIO/
  // in general or does PetscViewerBinaryOpen handle that?

  /* --- Output particle positions. --- */

  // Build the full file name.
  PetscCall(PetscStrcat(posfn, insert));
  PetscCall(PetscStrcat(posfn, ".bin"));

  // Create the binary viewer.
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, posfn, FILE_MODE_WRITE, &viewer));

  // Write particle positions to a binary file.
  PetscCall(DMSwarmCreateGlobalVectorFromField(swarmDM, DMSwarmPICField_coor, &target));
  PetscCall(PetscObjectSetName((PetscObject)target, "position"));
  PetscCall(VecView(target, viewer));
  PetscCall(DMSwarmDestroyGlobalVectorFromField(swarmDM, DMSwarmPICField_coor, &target));

  // Destroy the binary viewer.
  PetscCall(PetscViewerDestroy(&viewer));

  /* --- Output particle velocities. --- */

  // Build the full file name.
  PetscCall(PetscStrcat(velfn, insert));
  PetscCall(PetscStrcat(velfn, ".bin"));

  // Create the binary viewer.
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, velfn, FILE_MODE_WRITE, &viewer));

  // Write particle velocities to a binary file.
  PetscCall(DMSwarmCreateGlobalVectorFromField(swarmDM, "velocity", &target));
  PetscCall(PetscObjectSetName((PetscObject)target, "velocity"));
  PetscCall(VecView(target, viewer));
  PetscCall(DMSwarmDestroyGlobalVectorFromField(swarmDM, "velocity", &target));

  // Destroy the binary viewer.
  PetscCall(PetscViewerDestroy(&viewer));

  // Destroy the temporary vector object.
  PetscCall(VecDestroy(&target));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode OutputFluidHDF5(const char *insert, Context *ctx)
{
  char          name[PETSC_MAX_PATH_LEN]="fluid";
  PetscViewer   viewer;
  DM           *dms, dm;
  DM            fluidDM=ctx->fluidDM;
  PetscInt      Nf;
  char        **keys;
  PetscInt      field;
  Vec           target;
  Vec           moments=ctx->moments, rhs=ctx->potential.forcing, phi=ctx->potential.solution;

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Build the full file name.
  PetscCall(PetscStrcat(name, insert));
  PetscCall(PetscStrcat(name, ".hdf"));

  // Create the HDF5 viewer.
  PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &viewer));

  // Write fluid quantities to the HDF5 file.
  PetscCall(DMCreateFieldDecomposition(fluidDM, &Nf, &keys, NULL, &dms));
  for (field=0; field<Nf; field++) {
    dm = dms[field];
    PetscCall(DMGetGlobalVector(dm, &target));
    PetscCall(VecStrideGather(moments, field, target, INSERT_VALUES));
    PetscCall(PetscObjectSetName((PetscObject)target, keys[field]));
    PetscCall(VecView(target, viewer));
    PetscCall(DMRestoreGlobalVector(dm, &target));
  }

  // Release memory.
  for (field=0; field<Nf; field++) {
    PetscFree(keys[field]);
    PetscCall(DMDestroy(&dms[field]));
  }
  PetscFree(keys);
  PetscFree(dms);

  // Destroy the temporary vector object.
  PetscCall(VecDestroy(&target));

  // Write the forcing vector to the HDF5 file.
  PetscCall(PetscObjectSetName((PetscObject)rhs, "rhs"));
  PetscCall(VecView(rhs, viewer));

  // Write the electrostatic potential to the HDF5 file.
  PetscCall(PetscObjectSetName((PetscObject)phi, "potential"));
  PetscCall(VecView(phi, viewer));

  // Destroy the HDF5 viewer.
  PetscCall(PetscViewerDestroy(&viewer));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ViewLHS(KSP ksp, Context *ctx)
{
  PetscViewer viewer;
  Mat A, P;
  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  PetscCall(KSPGetOperators(ksp, &A, &P));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, "lhs.dat", FILE_MODE_WRITE, &viewer));
  PetscCall(MatView(A, viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}