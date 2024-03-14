# ECSPRT
Electrostatic Collisional Space Plasma Research Toolkit

This repository is still in the early development stage.

Until further notice, this repository requires building and compiling against
the `main` branch of PETSc in order to take advantage of crucial bug fixes.

ECSPRT (pronounced like "expert") currently comprises two main programs: a
plasma simulation and an electrostatic potential solver.

The ECSPRT plasma simulation is a hybrid particle-in-cell plasma code, meaning
it models electrons as a fluid and ions as particles. This allows it to avoid
certain temporal and spatial restrictions &mdash; namely, resolving the plasma
frequency and the Debye length &mdash; at the cost of failing to capture kinetic
electron physics. ECSPRT traces its lineage back to the EPPIC code developed
primarily at Boston University by Meers Oppenheim. It extends the 2D hybrid
capability of EPPIC to 3D by leveraging the Portable Extensible Toolkit for
Scientific Computation [(PETSc)](https://petsc.org/release/) for particle
communication (`DMSwarm`) and linear-system solvers (`KSP`).

The ECSPRT electrostatic potential solver provides a stand-alone application of
the same functionality that the ECSPRT plasma simulation uses to compute the
electrostatic potential at each time step. This stand-alone solver may be useful
for testing solution methods, as well as for computing the electrostatic
potential of an existing density distribution.

ECSPRT was created by Matt Young. This project welcomes contributions; it
especially encourages contributions from students, `git`/GitHub novices, and
software developers with limited knowledge of plasma physics.
