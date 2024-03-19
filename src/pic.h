#ifndef PIC_H
#define PIC_H

#include <petsc.h>

/* Defined types of density functions. */
typedef enum {
  DENSITY_FLAT_SOBOL,
  DENSITY_FLAT_REVERSE,
  DENSITY_FLAT_NORMAL,
  DENSITY_UNIFORM,
  DENSITY_UNIFORM_COORDINATES,
  DENSITY_UNIFORM_CENTERED,
  DENSITY_SINUSOIDAL,
  DENSITY_GAUSSIAN,
} DensityType;

/* Array that will hold names of density functions. */
extern const char *DensityTypes[];

/* Defined types of initial velocity distributions */
typedef enum {
  VELOCITIES_NORMAL
} VelocitiesType;

/* Array that will hold names of initial velocity distributions. */
extern const char *VelocitiesTypes[];

#endif // PIC_H
