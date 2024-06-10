# Changelog

## Developer Notes

When incrementing the version number to X.Y.Z, please do the following
* create a new subsection below **NEXT** with the title vX.Y.Z (YYYY-MM-DD)
* update the version number in `src/ecsprt.h`
* commit both files with the message "Increment version to X.Y.Z"
* create a tag named "vX.Y.Z" with the message "version X.Y.Z"
* push and follow tags via one of the following
  * `git push --follow-tags` on the command line
  * `Git: Push (Follow Tags)` from VS Code

## NEXT

## v0.6.0 (2024-06-10)

- Redefine boundary condition parameters
- Rename some CLI flags
- Add support for specifying number of particles per processor
- Deprecate use, and remove definition, of `NDIM`
- Fix memory bug and refactor to reduce memory leaks
- Allow user to pass options file as argument to executable
- Echo number of MPI processes at runtime
- Implement `--name` option for `build.sh`
- Add `make all` target
- Update `run.sh` diagnostic output
- Add color support to test script
- Make SLEPc features optional

## v0.5.0 (2024-03-21)

- Implement `-ndim {2,3}` option for `pic` simulations
- Significantly update and refactor initial particle distributions
- Extend tests

## v0.4.2 (2024-03-15)

- Reimplement particle rejection
- Fix bugs in `setup.sh` and `build.sh`
- Fix bugs in, and refactor, solver `main`
- Add initial tests

## v0.4.1 (2024-03-11)

- Fix bugs in particle initialization rejection function

## v0.4.0 (2024-03-07)

- Allow multiple `--options` arguments to `run.sh`
- Rename package to ECSPRT

## v0.3.1 (2024-03-07)

- Fix `--version` output bug
- Add text to `--version` output

## v0.3.0 (2024-03-07)

- Implement `--log-level` option for simulations
- Refactor runtime messages for simulations
- Make `-v/--verbose` more consistent in `run.sh`

## v0.2.1 (2024-03-07)

- Move DM viewing options to CLI only
- Fix bug in `build.sh` when user passes `--log none`
- Update runtime diagnostic messages

## v0.2.0 (2024-03-05)

- Fix bugs in and refactor `SetUpContext`
- Fix bug with multiple processors in `ComputeCollisions`
- Update coordinates in ions DM to physical values
- Significant modifications to `SobolDistribution`
- Add `--with-slepc-dir` option to `setup.sh`
- Various additional refactorings

## v0.1.1 (2023-12-12)

- Significant updates to `setup.sh`, `build.sh`, and `run.sh`
- Output final positions and velocities in binary format

## v0.1.0 (2023-12-08)

- Hello world!

