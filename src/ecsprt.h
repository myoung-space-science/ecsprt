#ifndef HYBRID_H
#define HYBRID_H

#include <petsc.h>
#include "constants.h"
#include "logging.h"
#include "global.h"

/* Version number, project full name, and project acronym.

These will appear when a user passes the --version flag.

See CHANGELOG.md for instructions on incrementing the version number.
*/
#define VERSION "0.6.0"
#define PROJECT "Electrostatic Collisional Space Plasma Research Toolkit"
#define ACRONYM "ECSPRT"

typedef enum {
  FORWARD,
  BACKWARD,
  CENTERED,
} DifferenceType;

extern const char *BCTypes[];
typedef enum {
  BC_PERIODIC,
  BC_INJECT_REFLECT,
  BC_INJECT_ADVECT,
} BCType;

/* MPI-specific parameters. */
typedef struct {
  PetscMPIInt rank; // global processor number
  PetscMPIInt size; // total number of processors
} MPI;

/* Parameters of the simulation spatial grid. */
typedef struct {
  PetscInt       Nx;  // global number of cells in x dimension
  PetscInt       Ny;  // global number of cells in y dimension
  PetscInt       Nz;  // global number of cells in z dimension
  PetscInt       nx;  // local number of cells in x dimension
  PetscInt       ny;  // local number of cells in y dimension
  PetscInt       nz;  // local number of cells in z dimension
  PetscReal      Lx;  // physical length of x dimension
  PetscReal      Ly;  // physical length of y dimension
  PetscReal      Lz;  // physical length of z dimension
  PetscReal      dx;  // physical cell spacing in x dimension
  PetscReal      dy;  // physical cell spacing in y dimension
  PetscReal      dz;  // physical cell spacing in z dimension
  PetscReal      x0;  // lower physical bound of x dimension
  PetscReal      y0;  // lower physical bound of y dimension
  PetscReal      z0;  // lower physical bound of z dimension
  PetscReal      x1;  // upper physical bound of x dimension
  PetscReal      y1;  // upper physical bound of y dimension
  PetscReal      z1;  // upper physical bound of z dimension
  DMBoundaryType xBC; // x-dimension boundary type
  DMBoundaryType yBC; // y-dimension boundary type
  DMBoundaryType zBC; // z-dimension boundary type
} Grid;

/* Bulk plasma parameters. */
typedef struct {
  PetscReal B0; // constant magnetic-field amplitude
  PetscReal E0; // constant electric-field amplitude
  PetscReal n0; // background density
} Plasma;

/* Parameters of a plasma species. */
typedef struct {
  PetscReal   q;       // charge
  PetscReal   m;       // mass
  PetscReal   nu;      // frequency of collisions with neutral particles
  PetscReal   Omega;   // gyrofrequency components
  PetscReal   kappa;   // magnetization components
  PetscReal   v0x;     // drift velocity in x dimension
  PetscReal   v0y;     // drift velocity in y dimension
  PetscReal   v0z;     // drift velocity in z dimension
  PetscReal   v0;      // drift-velocity magnitude
  PetscReal   vTx;     // thermal velocity in x dimension
  PetscReal   vTy;     // thermal velocity in y dimension
  PetscReal   vTz;     // thermal velocity in z dimension
  PetscReal   vT;      // thermal-velocity magnitude
  PetscReal   T;       // temperature
  BCType      xBC;     // x-axis boundary condition
  BCType      yBC;     // y-axis boundary condition
  BCType      zBC;     // z-axis boundary condition
  BCFunc      applyBC; // function to use when applying boundary conditions
  CollectFunc collect; // function to collect particle moments
  PushFunc    push;    // function to update particle velocities without collisions
  CollideFunc collide; // function to collide particles
} Species;

/* Parameters of the linear potential system. */
typedef struct {
  DM              dm;          // data manager for the logical grid
  Vec             solution;    // the solution vector
  Vec             forcing;     // the forcing vector
  LHSFunc         lhs;         // function to use when constructing LHS operator
  RHSFunc         rhs;         // function to use when constructing RHS function
  PetscReal       scale;       // value by which to scale RHS and LHS
  PetscBool       viewLHS;     // option to view LHS operator structure
  PetscInt        stencilSize; // number of points in the LHS-matrix stencil
  DMDAStencilType stencilType; // type of LHS-matrix stencil
} Linear;

/* Logging functions. */
typedef struct {
  LogFunction world;
  LogFunction self;
  LogFunction ranks;
  LogFunction checkpoint;
  LogFunction status;
} Loggers;

/* Parameters common to all applications. */
typedef struct {
  time_t    startTime;                      // global application start time
  time_t    endTime;                        // global application end time
  Loggers   log;                            // Logging functions
  MPI       mpi;                            // MPI parameters
  Grid      grid;                           // the spatial grid
  Species   electrons;                      // parameters of the electron fluid
  Species   ions;                           // parameters of the ion fluid
  Species   neutrals;                       // parameters of the neutral fluid
  Vec       moments;                        // fluid moments of the distribution function
  DM        fluidDM;                        // data manager for fluid quantities
  DM        swarmDM;                        // data manager for particles
  PetscReal gammaT;                         // electron thermal coefficient
  Plasma    plasma;                         // bulk plasma properties
  Linear    potential;                      // the electrostatic potential
  char      optionsLog[PETSC_MAX_PATH_LEN]; // path to output log
} Context;

/* Command-line options for all applications. */
typedef struct {
  PetscInt  ndim;    // number of spatial dimensions
  RHSType   rhsType; // type of RHS vector to use
  LHSType   lhsType; // type of LHS operator to use
  PetscInt  Nx;      // number of cells in x dimension
  PetscInt  Ny;      // number of cells in y dimension
  PetscInt  Nz;      // number of cells in z dimension
  PetscReal dx;      // physical cell spacing in x dimension
  PetscReal dy;      // physical cell spacing in y dimension
  PetscReal dz;      // physical cell spacing in z dimension
  PetscReal x0;      // lower physical bound of x dimension
  PetscReal y0;      // lower physical bound of y dimension
  PetscReal z0;      // lower physical bound of z dimension
  PetscReal x1;      // upper physical bound of x dimension
  PetscReal y1;      // upper physical bound of y dimension
  PetscReal z1;      // upper physical bound of z dimension
  PetscReal B0;      // constant magnetic-field amplitude
  PetscReal E0;      // constant electric-field amplitude
  PetscReal n0;      // plasma density
  PetscReal qi;      // ion charge
  PetscReal mi;      // ion mass
  PetscReal nui;     // frequency of collisions between ions and neutral particles
  PetscReal vi0x;    // ion drift velocity in x dimension
  PetscReal vi0y;    // ion drift velocity in y dimension
  PetscReal vi0z;    // ion drift velocity in z dimension
  PetscReal vi0;     // ion drift-velocity magnitude
  PetscReal viTx;    // ion thermal velocity in x dimension
  PetscReal viTy;    // ion thermal velocity in y dimension
  PetscReal viTz;    // ion thermal velocity in z dimension
  PetscReal viT;     // ion thermal-velocity magnitude
  PetscReal Ti;      // ion temperature
  PetscReal nue;     // frequency of collisions between electrons and neutral particles
  PetscReal ve0x;    // electron drift velocity in x dimension
  PetscReal ve0y;    // electron drift velocity in y dimension
  PetscReal ve0z;    // electron drift velocity in z dimension
  PetscReal ve0;     // electron drift-velocity magnitude
  PetscReal veTx;    // electron thermal velocity in x dimension
  PetscReal veTy;    // electron thermal velocity in y dimension
  PetscReal veTz;    // electron thermal velocity in z dimension
  PetscReal veT;     // electron thermal-velocity magnitude
  PetscReal Te;      // electron temperature
  PetscReal mn;      // neutral-particle mass
  PetscReal vn0x;    // neutral-particle drift velocity in x dimension
  PetscReal vn0y;    // neutral-particle drift velocity in y dimension
  PetscReal vn0z;    // neutral-particle drift velocity in z dimension
  PetscReal vn0;     // neutral-particle drift-velocity magnitude
  PetscReal vnTx;    // neutral-particle thermal velocity in x dimension
  PetscReal vnTy;    // neutral-particle thermal velocity in y dimension
  PetscReal vnTz;    // neutral-particle thermal velocity in z dimension
  PetscReal vnT;     // neutral-particle thermal-velocity magnitude
  PetscReal Tn;      // neutral-particle temperature
  BCType    xBC;     // x-axis boundary condition(s)
  BCType    yBC;     // y-axis boundary condition(s)
  BCType    zBC;     // z-axis boundary condition(s)
} CLI;

/* Fluid quantities defined at each spatial grid point. */
typedef struct {
  PetscReal n;
  PetscReal Gx;
  PetscReal Gy;
  PetscReal Gz;
} FluidNode;

#endif // HYBRID_H
