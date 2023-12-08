#ifndef CALCULUS_H
#define CALCULUS_H

#include <petsc.h>
#include "hybrid.h"

extern PetscErrorCode dFdx(PetscReal ***F, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type);
extern PetscErrorCode dFdy(PetscReal ***F, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type);
extern PetscErrorCode dFdz(PetscReal ***F, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type);

extern PetscErrorCode d2Fdxx(PetscReal ***F, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type);
extern PetscErrorCode d2Fdyy(PetscReal ***F, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type);
extern PetscErrorCode d2Fdzz(PetscReal ***F, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type);

extern PetscErrorCode d2Fdxy(PetscReal ***F, PetscReal dx, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType xtype, DifferenceType ytype);
extern PetscErrorCode d2Fdyx(PetscReal ***F, PetscReal dy, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ytype, DifferenceType xtype);
extern PetscErrorCode d2Fdxz(PetscReal ***F, PetscReal dx, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType xtype, DifferenceType ztype);
extern PetscErrorCode d2Fdzx(PetscReal ***F, PetscReal dz, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ztype, DifferenceType xtype);
extern PetscErrorCode d2Fdyz(PetscReal ***F, PetscReal dy, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ytype, DifferenceType ztype);
extern PetscErrorCode d2Fdzy(PetscReal ***F, PetscReal dz, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ztype, DifferenceType ytype);

extern PetscErrorCode DifferenceVector(PetscReal ***F, PetscReal x0, PetscReal y0, PetscReal z0, Grid grid, PetscReal f[NDIM]);
extern PetscErrorCode DotProduct(PetscReal a[NDIM], PetscReal b[NDIM], PetscReal *c);
extern PetscErrorCode CrossProduct(PetscReal a[NDIM], PetscReal b[NDIM], PetscReal c[NDIM]);

extern PetscErrorCode ComputeGradient(Vec F, Vec f[NDIM], Context *ctx);
extern PetscErrorCode ComputeLaplacian(Vec F, Vec *f, Context *ctx);

#endif // CALCULUS_H