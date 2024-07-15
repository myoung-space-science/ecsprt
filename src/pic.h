#ifndef PIC_H
#define PIC_H

#include <petsc.h>

/* Defined types of initial position distributions. */
typedef enum {
  PDIST_SOBOL,
  PDIST_LOCAL_SOBOL,
  PDIST_GLOBAL_SOBOL,
  PDIST_REVERSE,
  PDIST_NORMAL,
  PDIST_LOCAL_NORMAL,
  PDIST_GLOBAL_NORMAL,
  PDIST_UNIFORM,
  PDIST_SINUSOIDAL,
  PDIST_GAUSSIAN,
} PDistType;

/* Array that will hold names of density functions. */
extern const char *PDistTypes[];

/* Defined types of initial velocity distributions */
typedef enum {
  VDIST_NORMAL
} VDistType;

/* Array that will hold names of initial velocity distributions. */
extern const char *VDistTypes[];

#endif // PIC_H
