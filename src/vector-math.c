#include <petsc.h>
#include "hybrid.h"


/* Compute dF/dx to second order.

  Approximate f(i,j,k) = dF/dx at (i,j,k) using second-order centered, forward,
  or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dx [in]: The cell width. Passing dx <= 0.0 will cause this function to ignore
    this parameter and return the differential quantity dF. Doing so may be
    useful in cases when calling code prefers to scale the value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - type [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode dFdx(PetscReal ***F, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type)
{
  PetscReal v;

  PetscFunctionBeginUser;

  switch (type) {
    case CENTERED:
      v = F[k][j][i+1] - F[k][j][i-1];
      break;
    case FORWARD:
      v = -1.0*F[k][j][i+2] + 4.0*F[k][j][i+1] - 3.0*F[k][j][i];
      break;
    case BACKWARD:
      v = +3.0*F[k][j][i] - 4.0*F[k][j][i-1] + 1.0*F[k][j][i-2];
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }
  if (dx > 0.0) {
    *f = v / (2.0*dx);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute dF/dy to second order.

  Approximate f(i,j,k) = dF/dy at (i,j,k) using second-order centered, forward,
  or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dy [in]: The cell width. Passing dy <= 0.0 will cause this function to ignore
    this parameter and return the differential quantity dF. Doing so may be
    useful in cases when calling code prefers to scale the value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - type [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode dFdy(PetscReal ***F, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type)
{
  PetscReal v;

  PetscFunctionBeginUser;

  switch (type) {
    case CENTERED:
      v = F[k][j+1][i] - F[k][j-1][i];
      break;
    case FORWARD:
      v = -1.0*F[k][j+2][i] + 4.0*F[k][j+1][i] - 3.0*F[k][j][i];
      break;
    case BACKWARD:
      v = +3.0*F[k][j][i] - 4.0*F[k][j-1][i] + 1.0*F[k][j-2][i];
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }
  if (dy > 0.0) {
    *f = v / (2.0*dy);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute dF/dz to second order.

  Approximate f(i,j,k) = dF/dz at (i,j,k) using second-order centered, forward,
  or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dz [in]: The cell width. Passing dz <= 0.0 will cause this function to ignore
    this parameter and return the differential quantity dF. Doing so may be
    useful in cases when calling code prefers to scale the value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - type [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode dFdz(PetscReal ***F, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type)
{
  PetscReal v;

  PetscFunctionBeginUser;

  switch (type) {
    case CENTERED:
      v = F[k+1][j][i] - F[k-1][j][i];
      break;
    case FORWARD:
      v = -1.0*F[k+2][j][i] + 4.0*F[k+1][j][i] - 3.0*F[k][j][i];
      break;
    case BACKWARD:
      v = +3.0*F[k][j][i] - 4.0*F[k-1][j][i] + 1.0*F[k-2][j][i];
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }
  if (dz > 0.0) {
    *f = v / (2.0*dz);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dx)/dx to second order.

  Approximate f(i,j,k) = d(dF/dx)/dx at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dx [in]: The cell width. Passing dx <= 0.0 will cause this function to
    ignore this parameter and return the differential quantity dF. Doing so may
    be useful in cases when calling code prefers to scale the value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - type [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdxx(PetscReal ***F, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type)
{
  PetscReal v;

  PetscFunctionBeginUser;

  switch (type) {
    case CENTERED:
      v = F[k][j][i+1] - 2.0*F[k][j][i] + F[k][j][i-1];
      break;
    case FORWARD:
      v = -1.0*F[k][j][i+3] + 4.0*F[k][j][i+2] - 5.0*F[k][j][i+1] + 2.0*F[k][j][i];
      break;
    case BACKWARD:
      v = -1.0*F[k][j][i-3] + 4.0*F[k][j][i-2] - 5.0*F[k][j][i-1] + 2.0*F[k][j][i];
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }
  if (dx > 0.0) {
    *f = v / (dx*dx);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dy)/dy to second order.

  Approximate f(i,j,k) = d(dF/dy)/dy at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dy [in]: The cell width. Passing dy <= 0.0 will cause this function to
    ignore this parameter and return the differential quantity dF. Doing so may
    be useful in cases when calling code prefers to scale the value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - type [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdyy(PetscReal ***F, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type)
{
  PetscReal v;

  PetscFunctionBeginUser;

  switch (type) {
    case CENTERED:
      v = F[k][j+1][i] - 2.0*F[k][j][i] + F[k][j-1][i];
      break;
    case FORWARD:
      v = -1.0*F[k][j+3][i] + 4.0*F[k][j+2][i] - 5.0*F[k][j+1][i] + 2.0*F[k][j][i];
      break;
    case BACKWARD:
      v = -1.0*F[k][j-3][i] + 4.0*F[k][j-2][i] - 5.0*F[k][j-1][i] + 2.0*F[k][j][i];
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }
  if (dy > 0.0) {
    *f = v / (dy*dy);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dz)/dz to second order.

  Approximate f(i,j,k) = d(dF/dz)/dz at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dz [in]: The cell width. Passing dz <= 0.0 will cause this function to
    ignore this parameter and return the differential quantity dF. Doing so may
    be useful in cases when calling code prefers to scale the value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - type [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdzz(PetscReal ***F, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType type)
{
  PetscReal v;

  PetscFunctionBeginUser;

  switch (type) {
    case CENTERED:
      v = F[k+1][j][i] - 2.0*F[k][j][i] + F[k-1][j][i];
      break;
    case FORWARD:
      v = -1.0*F[k+3][j][i] + 4.0*F[k+2][j][i] - 5.0*F[k+1][j][i] + 2.0*F[k][j][i];
      break;
    case BACKWARD:
      v = -1.0*F[k-3][j][i] + 4.0*F[k-2][j][i] - 5.0*F[k-1][j][i] + 2.0*F[k][j][i];
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }
  if (dz > 0.0) {
    *f = v / (dz*dz);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dx)/dy to second order.

  Approximate f(i,j,k) = d(dF/dx)/dy at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dx, dy [in]: The cell widths. Passing dx <= 0.0 or dy <= 0.0 will cause this
    function to ignore this parameter and return the differential quantity dF.
    Doing so may be useful in cases when calling code prefers to scale the
    value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - xtype [in]: One of CENTERED, FOWARD, or BACKWARD.
  - ytype [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdxy(PetscReal ***F, PetscReal dx, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType xtype, DifferenceType ytype)
{
  PetscReal v;

  PetscFunctionBeginUser;

         if ((xtype == CENTERED) && (ytype == CENTERED)) {
    v = F[k][j+1][i+1] - F[k][j+1][i-1] - F[k][j-1][i+1] + F[k][j-1][i-1];
  } else if ((xtype == CENTERED) && (ytype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "CENTERED+FORWARD not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == CENTERED) && (ytype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "CENTERED+BACKWARD not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == FORWARD)  && (ytype == CENTERED)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+CENTERED not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == FORWARD)  && (ytype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+FORWARD not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == FORWARD)  && (ytype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+BACKWARD not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == BACKWARD) && (ytype == CENTERED)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+CENTERED not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == BACKWARD) && (ytype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+FORWARD not implemented for d(dF/dx)/dy.\n");
  } else if ((xtype == BACKWARD) && (ytype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+BACKWARD not implemented for d(dF/dx)/dy.\n");
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }

  if ((dx > 0.0) && (dy > 0.0)) {
    *f = v / (4.0*dx*dy);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dy)/dx to second order.

  Approximate f(i,j,k) = d(dF/dy)/dx at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dy, dx [in]: The cell widths. Passing dy <= 0.0 or dx <= 0.0 will cause this
    function to ignore this parameter and return the differential quantity dF.
    Doing so may be useful in cases when calling code prefers to scale the
    value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - ytype [in]: One of CENTERED, FOWARD, or BACKWARD.
  - xtype [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdyx(PetscReal ***F, PetscReal dy, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ytype, DifferenceType xtype)
{
  PetscFunctionBeginUser;
  PetscCall(d2Fdxy(F, dx, dy, i, j, k, f, xtype, ytype));
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dx)/dz to second order.

  Approximate f(i,j,k) = d(dF/dx)/dz at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dx, dz [in]: The cell widths. Passing dx <= 0.0 or dz <= 0.0 will cause this
    function to ignore this parameter and return the differential quantity dF.
    Doing so may be useful in cases when calling code prefers to scale the
    value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - xtype [in]: One of CENTERED, FOWARD, or BACKWARD.
  - ztype [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdxz(PetscReal ***F, PetscReal dx, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType xtype, DifferenceType ztype)
{
  PetscReal v;

  PetscFunctionBeginUser;

         if ((xtype == CENTERED) && (ztype == CENTERED)) {
    v = F[k+1][j][i+1] - F[k+1][j][i-1] - F[k-1][j][i+1] + F[k-1][j][i-1];
  } else if ((xtype == CENTERED) && (ztype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "CENTERED+FORWARD not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == CENTERED) && (ztype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "CENTERED+BACKWARD not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == FORWARD)  && (ztype == CENTERED)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+CENTERED not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == FORWARD)  && (ztype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+FORWARD not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == FORWARD)  && (ztype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+BACKWARD not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == BACKWARD) && (ztype == CENTERED)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+CENTERED not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == BACKWARD) && (ztype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+FORWARD not implemented for d(dF/dx)/dz.\n");
  } else if ((xtype == BACKWARD) && (ztype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+BACKWARD not implemented for d(dF/dx)/dz.\n");
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }

  if ((dx > 0.0) && (dz > 0.0)) {
    *f = v / (4.0*dx*dz);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dz)/dx to second order.

  Approximate f(i,j,k) = d(dF/dz)/dx at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dz, dx [in]: The cell widths. Passing dz <= 0.0 or dx <= 0.0 will cause this
    function to ignore this parameter and return the differential quantity dF.
    Doing so may be useful in cases when calling code prefers to scale the
    value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - ztype [in]: One of CENTERED, FOWARD, or BACKWARD.
  - xtype [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdzx(PetscReal ***F, PetscReal dz, PetscReal dx, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ztype, DifferenceType xtype)
{
  PetscFunctionBeginUser;
  PetscCall(d2Fdxz(F, dx, dz, i, j, k, f, xtype, ztype));
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dy)/dz to second order.

  Approximate f(i,j,k) = d(dF/dy)/dz at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dy, dz [in]: The cell widths. Passing dy <= 0.0 or dz <= 0.0 will cause this
    function to ignore this parameter and return the differential quantity dF.
    Doing so may be useful in cases when calling code prefers to scale the
    value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - ytype [in]: One of CENTERED, FOWARD, or BACKWARD.
  - ztype [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdyz(PetscReal ***F, PetscReal dy, PetscReal dz, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ytype, DifferenceType ztype)
{
  PetscReal v;

  PetscFunctionBeginUser;

         if ((ytype == CENTERED) && (ztype == CENTERED)) {
    v = F[k+1][j][i+1] - F[k+1][j][i-1] - F[k-1][j][i+1] + F[k-1][j][i-1];
  } else if ((ytype == CENTERED) && (ztype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "CENTERED+FORWARD not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == CENTERED) && (ztype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "CENTERED+BACKWARD not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == FORWARD)  && (ztype == CENTERED)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+CENTERED not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == FORWARD)  && (ztype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+FORWARD not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == FORWARD)  && (ztype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "FORWARD+BACKWARD not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == BACKWARD) && (ztype == CENTERED)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+CENTERED not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == BACKWARD) && (ztype == FORWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+FORWARD not implemented for d(dF/dy)/dz.\n");
  } else if ((ytype == BACKWARD) && (ztype == BACKWARD)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "BACKWARD+BACKWARD not implemented for d(dF/dy)/dz.\n");
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Accepted difference types are FORWARD, BACKWARD, and CENTERED.\n");
  }

  if ((dy > 0.0) && (dz > 0.0)) {
    *f = v / (4.0*dy*dz);
  } else {
    *f = v;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute d(dF/dz)/dy to second order.

  Approximate f(i,j,k) = d(dF/dz)/dy at (i,j,k) using second-order centered,
  forward, or backward differencing scheme.

  - F [in]: The antiderivative function.
  - dz, dy [in]: The cell widths. Passing dz <= 0.0 or dy <= 0.0 will cause this
    function to ignore this parameter and return the differential quantity dF.
    Doing so may be useful in cases when calling code prefers to scale the
    value.
  - i [in]: The x-axis index.
  - j [in]: The y-axis index.
  - k [in]: The z-axis index.
  - f [out]: A pointer to the scalar result.
  - ztype [in]: One of CENTERED, FOWARD, or BACKWARD.
  - ytype [in]: One of CENTERED, FOWARD, or BACKWARD.

*/
PetscErrorCode d2Fdzy(PetscReal ***F, PetscReal dz, PetscReal dy, PetscInt i, PetscInt j, PetscInt k, PetscReal *f, DifferenceType ztype, DifferenceType ytype)
{
  PetscFunctionBeginUser;
  PetscCall(d2Fdyz(F, dy, dz, i, j, k, f, ytype, ztype));
  PetscFunctionReturn(PETSC_SUCCESS);
}


// TODO: Refactor this with difference functions above.
/* Compute a vector of central differences from F.

This function was designed as the first step in computing the gradient of the
scalar function F(x, y, z) at (x0, y0, z0). It computes the numerator of each
finite-difference term using 2nd-order centered, forward, and backward
approximations. It assumes that F contains the appropriate ghost nodes.
*/
PetscErrorCode DifferenceVector(PetscReal ***F, PetscReal x0, PetscReal y0, PetscReal z0, Grid grid, PetscReal f[NDIM])
{
  PetscInt    Nx=grid.Nx, Ny=grid.Ny, Nz=grid.Nz;
  PetscInt    ixl, ixh, iyl, iyh, izl, izh;
  PetscReal   wxh, wyh, wzh;
  PetscReal   hhh, lhh, hlh, llh, hhl, lhl, hll, lll;
  PetscReal   whh, whl, wlh, wll, Ewh, Ewl;

  PetscFunctionBeginUser;

  // Compute the x-dimension neighbors and corresponding weights.
  ixl = (PetscInt)x0;
  ixh = ixl+1;
  wxh = x0 - (PetscReal)ixl;
  // Compute the y-dimension neighbors and corresponding weights.
  iyl = (PetscInt)y0;
  iyh = iyl+1;
  wyh = y0 - (PetscReal)iyl;
  // Compute the z-dimension neighbors and corresponding weights.
  izl = (PetscInt)z0;
  izh = izl+1;
  wzh = z0 - (PetscReal)izl;
  // Compute the central difference in x at each grid point.
  if (ixl >= 0) {
    // 2nd-order central difference at ixl
    hhl = F[izh][iyh][ixl+1] - F[izh][iyh][ixl-1];
    lhl = F[izl][iyh][ixl+1] - F[izl][iyh][ixl-1];
    hll = F[izh][iyl][ixl+1] - F[izh][iyl][ixl-1];
    lll = F[izl][iyl][ixl+1] - F[izl][iyl][ixl-1];
  } else {
    // 2nd-order forward difference at ixl
    hhl = -1.0*F[izh][iyh][ixl+2] + 4.0*F[izh][iyh][ixl+1] - 3.0*F[izh][iyh][ixl];
    lhl = -1.0*F[izl][iyh][ixl+2] + 4.0*F[izl][iyh][ixl+1] - 3.0*F[izl][iyh][ixl];
    hll = -1.0*F[izh][iyl][ixl+2] + 4.0*F[izh][iyl][ixl+1] - 3.0*F[izh][iyl][ixl];
    lll = -1.0*F[izl][iyl][ixl+2] + 4.0*F[izl][iyl][ixl+1] - 3.0*F[izl][iyl][ixl];
  }
  if (ixh < Nx) {
    // 2nd-order central difference at ixh
    hhh = F[izh][iyh][ixh+1] - F[izh][iyh][ixh-1];
    lhh = F[izl][iyh][ixh+1] - F[izl][iyh][ixh-1];
    hlh = F[izh][iyl][ixh+1] - F[izh][iyl][ixh-1];
    llh = F[izl][iyl][ixh+1] - F[izl][iyl][ixh-1];
  } else {
    // 2nd-order backward difference at ixh
    hhh = +3.0*F[izh][iyh][ixh] - 4.0*F[izh][iyh][ixh-1] + 1.0*F[izh][iyh][ixh-2];
    lhh = +3.0*F[izl][iyh][ixh] - 4.0*F[izl][iyh][ixh-1] + 1.0*F[izl][iyh][ixh-2];
    hlh = +3.0*F[izh][iyl][ixh] - 4.0*F[izh][iyl][ixh-1] + 1.0*F[izh][iyl][ixh-2];
    llh = +3.0*F[izl][iyl][ixh] - 4.0*F[izl][iyl][ixh-1] + 1.0*F[izl][iyl][ixh-2];
  }
  whh = hlh + wyh*(hhh - hlh);
  whl = hll + wyh*(hhl - hll);
  wlh = llh + wyh*(lhh - llh);
  wll = lll + wyh*(lhl - lll);
  Ewh = wlh + wxh*(whh - wlh);
  Ewl = wll + wxh*(whl - wll);
  f[0] = Ewl + wzh*(Ewh - Ewl);
  // Compute the central difference in y at each grid point.
  if (iyl >= 0) {
    // 2nd-order central difference at iyl
    hlh = F[izh][iyl+1][ixh] - F[izh][iyl-1][ixh];
    llh = F[izl][iyl+1][ixh] - F[izl][iyl-1][ixh];
    hll = F[izh][iyl+1][ixl] - F[izh][iyl-1][ixl];
    lll = F[izl][iyl+1][ixl] - F[izl][iyl-1][ixl];
  } else {
    // 2nd-order forward difference at iyl
    hlh = -1.0*F[izh][iyl+2][ixh] + 4.0*F[izh][iyl+1][ixh] - 3.0*F[izh][iyl][ixh];
    llh = -1.0*F[izl][iyl+2][ixh] + 4.0*F[izl][iyl+1][ixh] - 3.0*F[izl][iyl][ixh];
    hll = -1.0*F[izh][iyl+2][ixl] + 4.0*F[izh][iyl+1][ixl] - 3.0*F[izh][iyl][ixl];
    lll = -1.0*F[izl][iyl+2][ixl] + 4.0*F[izl][iyl+1][ixl] - 3.0*F[izl][iyl][ixl];
  }
  if (iyh < Ny) {
    // 2nd-order central difference at iyh
    hhh = F[izh][iyh+1][ixh] - F[izh][iyh-1][ixh];
    lhh = F[izl][iyh+1][ixh] - F[izl][iyh-1][ixh];
    hhl = F[izh][iyh+1][ixl] - F[izh][iyh-1][ixl];
    lhl = F[izl][iyh+1][ixl] - F[izl][iyh-1][ixl];
  } else {
    // 2nd-order backward difference at iyh
    hhh = +3.0*F[izh][iyh][ixh] - 4.0*F[izh][iyh-1][ixh] + 1.0*F[izh][iyh-2][ixh];
    lhh = +3.0*F[izl][iyh][ixh] - 4.0*F[izl][iyh-1][ixh] + 1.0*F[izl][iyh-2][ixh];
    hhl = +3.0*F[izh][iyh][ixl] - 4.0*F[izh][iyh-1][ixl] + 1.0*F[izh][iyh-2][ixl];
    lhl = +3.0*F[izl][iyh][ixl] - 4.0*F[izl][iyh-1][ixl] + 1.0*F[izl][iyh-2][ixl];
  }
  whh = hlh + wyh*(hhh - hlh);
  whl = hll + wyh*(hhl - hll);
  wlh = llh + wyh*(lhh - llh);
  wll = lll + wyh*(lhl - lll);
  Ewh = wlh + wxh*(whh - wlh);
  Ewl = wll + wxh*(whl - wll);
  f[1] = Ewl + wzh*(Ewh - Ewl);
  // Compute the central difference in z at each grid point.
  if (izl >= 0) {
    // 2nd-order central difference at izl
    lhh = F[izl+1][iyh][ixh] - F[izl-1][iyh][ixh];
    llh = F[izl+1][iyl][ixh] - F[izl-1][iyl][ixh];
    lhl = F[izl+1][iyh][ixl] - F[izl-1][iyh][ixl];
    lll = F[izl+1][iyl][ixl] - F[izl-1][iyl][ixl];
  } else {
    // 2nd-order forward difference at izl
    lhh = -1.0*F[izl+2][iyh][ixh] + 4.0*F[izl+1][iyh][ixh] - 3.0*F[izl][iyh][ixh];
    llh = -1.0*F[izl+2][iyl][ixh] + 4.0*F[izl+1][iyl][ixh] - 3.0*F[izl][iyl][ixh];
    lhl = -1.0*F[izl+2][iyh][ixl] + 4.0*F[izl+1][iyh][ixl] - 3.0*F[izl][iyh][ixl];
    lll = -1.0*F[izl+2][iyl][ixl] + 4.0*F[izl+1][iyl][ixl] - 3.0*F[izl][iyl][ixl];
  }
  if (izh < Nz) {
    // 2nd-order central difference at izh
    hhh = F[izh+1][iyh][ixh] - F[izh-1][iyh][ixh];
    hlh = F[izh+1][iyl][ixh] - F[izh-1][iyl][ixh];
    hhl = F[izh+1][iyh][ixl] - F[izh-1][iyh][ixl];
    hll = F[izh+1][iyl][ixl] - F[izh-1][iyl][ixl];
  } else {
    // 2nd-order backward difference at izh
    hhh = +3.0*F[izh][iyh][ixh] - 4.0*F[izh-1][iyh][ixh] + 1.0*F[izh-2][iyh][ixh];
    hlh = +3.0*F[izh][iyl][ixh] - 4.0*F[izh-1][iyl][ixh] + 1.0*F[izh-2][iyl][ixh];
    hhl = +3.0*F[izh][iyh][ixl] - 4.0*F[izh-1][iyh][ixl] + 1.0*F[izh-2][iyh][ixl];
    hll = +3.0*F[izh][iyl][ixl] - 4.0*F[izh-1][iyl][ixl] + 1.0*F[izh-2][iyl][ixl];
  }
  whh = hlh + wyh*(hhh - hlh);
  whl = hll + wyh*(hhl - hll);
  wlh = llh + wyh*(lhh - llh);
  wll = lll + wyh*(lhl - lll);
  Ewh = wlh + wxh*(whh - wlh);
  Ewl = wll + wxh*(whl - wll);
  f[2] = Ewl + wzh*(Ewh - Ewl);

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the dot product given by $\vec{c} = \vec{a} \cdot \vec{b}$. */
PetscErrorCode DotProduct(PetscReal a[NDIM], PetscReal b[NDIM], PetscReal *c)
{
  PetscFunctionBeginUser;

  *c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

  PetscFunctionReturn(PETSC_SUCCESS);
}


/* Compute the cross product given by $\vec{c} = \vec{a} \times \vec{b}$. */
PetscErrorCode CrossProduct(PetscReal a[NDIM], PetscReal b[NDIM], PetscReal c[NDIM])
{
  PetscFunctionBeginUser;

  c[0] = a[1]*b[2] - b[2]*a[1];
  c[1] = a[2]*b[0] - b[0]*a[2];
  c[2] = a[0]*b[1] - b[1]*a[0];

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeGradient(Vec F, Vec f[NDIM], Context *ctx)
{
  PetscReal    dx=ctx->grid.dx;
  PetscReal    dy=ctx->grid.dy;
  PetscReal    dz=ctx->grid.dz;
  DM           fluidDM=ctx->fluidDM;
  PetscInt     i0, j0, k0;
  PetscInt     ni, nj, nk;
  PetscInt     i, j, k;
  PetscReal    hx, hy, hz;
  Mat          Ax, Ay, Az;
  MatStencil   row;
  PetscInt     stencilSize=2;
  PetscReal    vals[stencilSize];
  MatStencil   cols[stencilSize];

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Compute differential scale factors.
  hx = 1.0 / (2.0*dx);
  hy = 1.0 / (2.0*dy);
  hz = 1.0 / (2.0*dz);

  // Create the operator matrices.
  PetscCall(DMCreateMatrix(fluidDM, &Ax));
  PetscCall(DMCreateMatrix(fluidDM, &Ay));
  PetscCall(DMCreateMatrix(fluidDM, &Az));

  // Compute the discrete derivative in each dimension.
  // Get this processor's indices.
  PetscCall(DMDAGetCorners(fluidDM, &i0, &j0, &k0, &ni, &nj, &nk));

  // Loop over grid points. [DEV] Assume periodic BC.
  for (k=k0; k<k0+nk; k++) {
    for (j=j0; j<j0+nj; j++) {
      for (i=i0; i<i0+ni; i++) {
        row.i = i; row.j = j; row.k = k;
        vals[0] = -hx;
        cols[0].i = i-1;
        cols[0].j = j;
        cols[0].k = k;
        vals[1] = +hx;
        cols[1].i = i+1;
        cols[1].j = j;
        cols[1].k = k;
        PetscCall(MatSetValuesStencil(Ax, 1, &row, stencilSize, cols, vals, INSERT_VALUES));
        vals[0] = -hy;
        cols[0].i = i;
        cols[0].j = j-1;
        cols[0].k = k;
        vals[1] = +hy;
        cols[1].i = i;
        cols[1].j = j+1;
        cols[1].k = k;
        PetscCall(MatSetValuesStencil(Ay, 1, &row, stencilSize, cols, vals, INSERT_VALUES));
        vals[0] = -hz;
        cols[0].i = i;
        cols[0].j = j;
        cols[0].k = k-1;
        vals[1] = +hz;
        cols[1].i = i;
        cols[1].j = j;
        cols[1].k = k+1;
        PetscCall(MatSetValuesStencil(Az, 1, &row, stencilSize, cols, vals, INSERT_VALUES));
      }
    }
  }

  // Assemble the distributed operator matrices.
  PetscCall(MatAssemblyBegin(Ax, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(Ay, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(Az, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Ax, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Ay, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Az, MAT_FINAL_ASSEMBLY));

  // Compute components of the discrete gradient.
  PetscCall(MatMult(Ax, F, f[0]));
  PetscCall(MatMult(Ay, F, f[1]));
  PetscCall(MatMult(Az, F, f[2]));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeLaplacian(Vec F, Vec *f, Context *ctx)
{
  PetscReal    dx=ctx->grid.dx;
  PetscReal    dy=ctx->grid.dy;
  PetscReal    dz=ctx->grid.dz;
  DM           fluidDM=ctx->fluidDM;
  PetscInt     i0, j0, k0;
  PetscInt     ni, nj, nk;
  PetscInt     i, j, k;
  PetscReal    hxx, hyy, hzz;
  PetscReal    vpjk, vmjk, vipk, vimk, vijp, vijm;
  PetscReal    vijk;
  Mat          A;
  MatStencil   row;
  PetscInt     stencilSize=7;
  PetscReal    vals[stencilSize];
  MatStencil   cols[stencilSize];

  PetscFunctionBeginUser;
  ctx->log.checkpoint("\n--> Entering %s <--\n", __func__);

  // Compute differential scale factors.
  hxx =  1.0 / (dx*dx);
  hyy =  1.0 / (dy*dy);
  hzz =  1.0 / (dz*dz);

  // Assign the star-stencil coefficients.
  vpjk =  hxx;
  vmjk =  hxx;
  vipk =  hyy;
  vimk =  hyy;
  vijp =  hzz;
  vijm =  hzz;

  // Assign the diagonal coefficient.
  vijk = -(vpjk + vipk + vijp + vmjk + vimk + vijm);

  // Create the operator matrix.
  PetscCall(DMCreateMatrix(fluidDM, &A));

  // Get this processor's indices.
  PetscCall(DMDAGetCorners(fluidDM, &i0, &j0, &k0, &ni, &nj, &nk));

  // Loop over grid points. [DEV] Assume periodic BC.
  for (k=k0; k<k0+nk; k++) {
    for (j=j0; j<j0+nj; j++) {
      for (i=i0; i<i0+ni; i++) {
        row.i = i; row.j = j; row.k = k;
        // Assign the value at node (i+1, j, k)
        vals[0] = vpjk;
        cols[0].i = i+1;
        cols[0].j = j;
        cols[0].k = k;
        // Assign the value at node (i-1, j, k)
        vals[1] = vmjk;
        cols[1].i = i-1;
        cols[1].j = j;
        cols[1].k = k;
        // Assign the value at node (i, j+1, k)
        vals[2] = vipk;
        cols[2].i = i;
        cols[2].j = j+1;
        cols[2].k = k;
        // Assign the value at node (i, j-1, k)
        vals[3] = vimk;
        cols[3].i = i;
        cols[3].j = j-1;
        cols[3].k = k;
        // Assign the value at node (i, j, k+1)
        vals[4] = vijp;
        cols[4].i = i;
        cols[4].j = j;
        cols[4].k = k+1;
        // Assign the value at node (i, j, k-1)
        vals[5] = vijm;
        cols[5].i = i;
        cols[5].j = j;
        cols[5].k = k-1;
        // Assign the value at node (i, j, k)
        vals[6] = vijk;
        cols[6].i = i;
        cols[6].j = j;
        cols[6].k = k;
        PetscCall(MatSetValuesStencil(A, 1, &row, stencilSize, cols, vals, INSERT_VALUES));
      }
    }
  }

  // Assemble the distributed operator matrix.
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  // Compute the discrete Laplacian.
  PetscCall(MatMult(A, F, *f));

  ctx->log.checkpoint("\n--> Exiting %s <--\n\n", __func__);
  PetscFunctionReturn(PETSC_SUCCESS);
}


