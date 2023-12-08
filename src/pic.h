#ifndef SIMULATION_H
#define SIMULATION_H

#include <petsc.h>

/* Supported density functions. */
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

extern const char *DensityTypes[];

#endif // SIMULATION_H
