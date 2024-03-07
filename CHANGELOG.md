# Changelog

## Developer Notes

When incrementing the version number to X.Y.Z, please do the following
* create a new subsection here (below **NEXT**) with the title vX.Y.Z (YYYY-MM-DD)
* update the version number in `src/hybrid.h`
* create a tag named "vX.Y.Z" with the message "version X.Y.Z"

## NEXT

## v0.2.1

- Move DM viewing options to CLI only
- Fix bug in `build.sh` when user passes `--log none`
- Update runtime diagnostic messages

## v0.2.0

- Fix bugs in and refactor `SetUpContext`
- Fix bug with multiple processors in `ComputeCollisions`
- Update coordinates in ions DM to physical values
- Significant modifications to `SobolDistribution`
- Add `--with-slepc-dir` option to `setup.sh`
- Various additional refactorings

## v0.1.1

- Significant updates to `setup.sh`, `build.sh`, and `run.sh`
- Output final positions and velocities in binary format

## v0.1.0

- Hello world!

