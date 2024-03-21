#!/bin/bash

cwd=$(dirname $(realpath $0))
rundir=${cwd}/runs

name=pic-gmres-hypre-3D
for n in 1 2; do
    ./run.sh -n $n -p pic -r $rundir -o $name \
        --options ${cwd}/grid3D.ini \
        --options ${cwd}/gmres-hypre.ini \
        --options ${cwd}/fb-physics.ini \
        --require-options \
        args \
        -Nt 2 \
        -dt 1e-5
    success=$?
    # TODO: Further analyze output before declaring success.
    if [ $success == 0 ]; then
        echo "Test $name passed on $n processors"
        rm -rf $rundir/$name
        rm -f $rundir/latest
    else
        echo "Test $name failed on $n processors"
        exit 1
    fi
done

name=solver-gmres-hypre-3D
for n in 1 2; do
    ./run.sh -n $n -p solver -r $rundir -o $name \
        --options ${cwd}/grid3D.ini \
        --options ${cwd}/gmres-hypre.ini \
        --options ${cwd}/fb-physics.ini \
        --require-options \
        args \
        --input ~/ecsprt/test-data/cos_xyz_32x32x32.h5 \
        --lhs-eigenvalues \
        --flux-scale 100.0
    success=$?
    # TODO: Further analyze output before declaring success.
    if [ $success == 0 ]; then
        echo "Test $name passed on $n processors"
        rm -rf $rundir/$name
        rm -f $rundir/latest
    else
        echo "Test $name failed on $n processors"
        exit 1
    fi
done

# TODO:
# - test solver very small case with --view-lhs
# - test solver small case with --lhs-eigenvalues
# - test all cases in 2D when possible
