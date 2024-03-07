#!/bin/bash

# See https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
# - e => Exit immediately if a pipeline returns non-zero status.
# - u => Treat unset variables and parameters (with certain exceptions) as
#        errors when performing parameter expansion.
set -eu

# Define special flag to indicate success.
success="@!SUCCESS!@"

# Define special unicode characters.
u_warning="\U1F937"
u_success="\u2728"
u_failure="\U1F62d"

# Declare the name of the target program.
prog=
cli=${0##*/}

# Declare build-related directories and files.
cwd=$(pwd)
srcdir=${cwd}/src
bindir=${cwd}/bin
mkdir -p "${bindir}"
buildlog=${cwd}/build.log
> ${buildlog}

# Set option defaults.
verbose=0
debug=0
optimize=0
from_clean=0
petsc_dir=
slepc_dir=
petsc_arch=

# Define text formatting commands.
# - textbf: Use bold-face.
# - textnm: Reset to normal.
# - startul: Start underlined text.
# - endul: End underlined text.
textbf=$(tput bold)
textnm=$(tput sgr0)
startul=$(tput smul)
endul=$(tput rmul)

# This is the CLI's main help text.
show_help()
{
    echo "
${textbf}NAME${textnm}
        $cli - Build the hybrid simulation or potential solver

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli${textnm} -p {pic,solver} [${startul}OPTIONS${endul}]

${textbf}DESCRIPTION${textnm}
        This script will perform all the steps to build TARGET.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-p${textnm}, ${textbf}--program${textnm}=${startul}PROG${endul}
                The target program (required).
        ${textbf}--from-clean${textnm}
                Run make clean in the target src directory before building the executable.
        ${textbf}--debug${textnm}
                Compile with debugging flags.
        ${textbf}--optimize${textnm}
                Compile with optimization flags.
        ${textbf}--with-petsc-dir DIR${textnm}
                Use DIR as PETSC_DIR when linking to PETSc.
        ${textbf}--with-slepc-dir DIR${textnm}
                Use DIR as SLEPC_DIR when linking to SLEPc.
        ${textbf}--with-petsc-arch ARCH${textnm}
                Use ARCH as PETSC_ARCH when linking to PETSc and SLEPc.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print informative messages during the build process.

${textbf}NOTES${textnm}
        This script will use values of PETSC_DIR, SLEPC_DIR, and PETSC_ARCH 
        set by the enviroment in the absence of the corresponding --with* options.
        If passed, the given value will override the environment value, if any.

"
}

# This will run if the CLI gets an unrecognized option.
report_bad_arg()
{
    printf "\nUnrecognized command: ${1}\n\n"
}

# Read command-line arguments.
TEMP=$(getopt \
    -n 'build.sh' \
    -o 'hvp:' \
    -l 'help,verbose' \
    -l 'program:' \
    -l 'from-clean' \
    -l 'debug,optimize' \
    -l 'with-petsc-dir:,with-slepc-dir:,with-petsc-arch:' \
    -- "$@")

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to parse command-line arguments."
    exit 1
fi

eval set -- "$TEMP"
unset TEMP

while [ $# -gt 0 ]; do
    case "$1" in
        '-h'|'--help')
            show_help
            exit
        ;;
        '-p'|'--program')
            prog="${2}"
            shift 2
            continue
        ;;
        '--from-clean')
            from_clean=1
            shift
            continue
        ;;
        '--debug')
            debug=1
            shift
            continue
        ;;
        '--optimize')
            optimize=1
            shift
            continue
        ;;
        '--with-petsc-dir')
            petsc_dir="${2}"
            shift 2
            continue
        ;;
        '--with-slepc-dir')
            slepc_dir="${2}"
            shift 2
            continue
        ;;
        '--with-petsc-arch')
            petsc_arch="${2}"
            shift 2
            continue
        ;;
        '-v'|'--verbose')
            verbose=1
            shift
            continue
        ;;
        '--')
            shift
            break
        ;;
    esac
done

# Alert the user to unrecognized arguments.
extra_args="$@"
if [ -n "${extra_args}" ]; then
    show_help
    echo
    echo -e "${u_warning} Got the following unrecognized argument(s):"
    echo
    for arg in ${extra_args[@]}; do
        echo ">   $arg"
    done
    echo
    echo "Possible causes of this error:"
    echo "- This script does not accept any positional arguments. "
    echo "  Please pass the required program name via the -p/--program option."
    echo "- This script does not support passing arbitrary arguments to ${textbf}make${textnm}."
    echo "  If you want (and know how to) build the code in a way that this script does not provide, "
    echo "  you must directly call ${textbf}make${textnm} with appropriate arguments."
    echo
    exit 1
fi

stage=

cleanup() {
    if [ "${stage}" != "${success}" ]; then
        if [ -n "${stage}" ]; then
            exitmsg="${stage} stage failed. See ${buildlog} for details."
        else
            exitmsg="Unknown error."
        fi
        echo 
        echo -e ${u_failure} $exitmsg
        exit 1
    else
        echo 
        echo -e "${u_success} Success! ${u_success}"
        echo 
        if [ $debug == 1 ]; then
            echo "For debugging on a single processor, try"
            echo "$ gdb --args ${bindir}/${prog} [ARGS]"
        else
            echo "To run on a single processor, try"
            echo "$ ${bindir}/${prog} [ARGS]"
            echo 
            echo "To run on N processors, try"
            echo "$ mpiexec -n N ${bindir}/${prog} [ARGS]"
        fi
    fi
}

trap cleanup EXIT

mark_stage() {
    stage="${1}"
    if [ $verbose == 1 ]; then
        echo "[${cli}] Executing ${stage} stage"
    fi
}

print_width() {
    local c=${1}

    printf "=%.0s" $(seq 1 $COLUMNS)
    echo 
}

# Set PETSc and SLEPc paths to user values, if given.
if [ -n "${petsc_dir}" ]; then
    export PETSC_DIR=${petsc_dir}
fi
if [ -n "${slepc_dir}" ]; then
    export SLEPC_DIR=${slepc_dir}
fi
if [ -n "${petsc_arch}" ]; then
    export PETSC_ARCH=${petsc_arch}
fi

stage="setup"

# Refuse to run without a target program.
if [ -z "${prog}" ]; then
    show_help &>> ${buildlog}
    echo
    echo "ERROR: Missing target program." &>> ${buildlog}
    exit 1
fi

if [ ${from_clean} == 1 ]; then
    # Mark this stage.
    mark_stage "clean"
    pushd "${srcdir}" &> /dev/null \
    && make clean &>> ${buildlog} \
    && popd &> /dev/null
fi

# Mark this stage.
mark_stage "build"

# Raise an error for empty PETSc or SLEPc path variables.
petsc_slepc_vars=(
    PETSC_DIR
    SLEPC_DIR
    PETSC_ARCH
)
for var in ${petsc_slepc_vars[@]}; do
    if [ -v ${var} ]; then
        echo "Using ${var}=${!var}" &>> ${buildlog}
    else
        echo "ERROR: ${var} is undefined." &>> ${buildlog}
        exit 1
    fi
done

# Build the executable in the source directory.
pushd ${srcdir} &> /dev/null
cppflags=
cflags=
progext=
if [ ${debug} == 1 ]; then
    cppflags="${cppflags} -DDEBUG"
    cflags="${cflags} -g -O0"
    progext="-dbg"
fi
if [ ${optimize} == 1 ]; then
    cppflags="${cppflags} -DNDEBUG"
    cflags="${cflags} -O3"
    progext="-opt"
fi
echo &>> ${buildlog}
print_width "=" &>> ${buildlog}
echo "    Running make on ECSPERT" &>> ${buildlog}
print_width "=" &>> ${buildlog}
make ${prog} CPPFLAGS="${cppflags}" CFLAGS="${cflags}" &>> ${buildlog}
popd &> /dev/null

if [ -n "${progext}" ]; then
    pushd ${bindir} &> /dev/null \
    && /bin/cp ${prog} ${prog}${progext} \
    && popd &> /dev/null
fi

# Signal success to the clean-up function.
stage=$success

