#!/bin/bash

cwd=$(dirname $(realpath $0))
rundir=${cwd}/runs

# TODO:
# - test solver small case with --lhs-eigenvalues
# - test solver very small case with --view-lhs

for N in 2 3; do
    for n in 1 2; do
        for dist in sobol sinusoidal; do

            name="pic-gmres-hypre-${N}D-$dist"
            ./run.sh -n $n -p pic -r $rundir -o $name \
                --options ${cwd}/grid3D.ini \
                --options ${cwd}/gmres-hypre.ini \
                --options ${cwd}/fb-physics.ini \
                --require-options \
                args \
                -ndim $N \
                -Nt 2 \
                -dt 1e-5
            success=$?
            # TODO: Further analyze output before declaring success.
            if [ $success == 0 ]; then
                echo "Test $name passed on $n processor(s)"
                rm -rf $rundir/$name
                rm -f $rundir/latest
            else
                echo "Test $name failed on $n processor(s)"
                exit 1
            fi

            name="solver-gmres-hypre-${N}D-$dist"
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
                echo "Test $name passed on $n processor(s)"
                rm -rf $rundir/$name
                rm -f $rundir/latest
            else
                echo "Test $name failed on $n processor(s)"
                exit 1
            fi

        done
    done
done

