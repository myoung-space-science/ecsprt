# ECSPERT
The Electrostatic Collisional Simulation of Plasma Evolution, Redistribution, and Transport

This repository is still in the early development stage.

ECSPERT is a hybrid particle-in-cell plasma code, meaning it models electrons as
a fluid and ions as particles. This allows it to avoid certain temporal and
spatial restrictions &mdash; namely, resolving the plasma frequency and the
Debye length &mdash; at the cost of failing to capture kinetic electron physics.
ECSPERT traces its lineage back to the EPPIC code developed primarily at Boston
University by Meers Oppenheim. It extends the 2D hybrid capability of EPPIC to
3D by leveraging the Portable Extensible Toolkit for Scientific Computation
[(PETSc)](https://petsc.org/release/) for particle communication (`DMSwarm`) and
linear-system solvers (`KSP`).

ECSPERT was created by Matt Young. This project welcomes contributions; it
especially encourages contributions from students, `git`/GitHub novices, and
software developers with limited knowledge of plasma physics.
