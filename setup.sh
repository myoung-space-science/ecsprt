# See https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
# - e => Exit immediately if a pipeline returns non-zero status.
# - u => Treat unset variables and parameters (with certain exceptions) as
#   errors when performing parameter expansion.
set -eu

# Define special flag to indicate success.
success="@!SUCCESS!@"

# Define special unicode characters.
u_warning="\U1F937"
u_success="\u2728"
u_error="\U1F62d"

# Define text formatting commands.
# - textbf: Use bold-face.
# - textnm: Reset to normal.
# - startul: Start underlined text.
# - endul: End underlined text.
textbf=$(tput bold)
textnm=$(tput sgr0)
startul=$(tput smul)
endul=$(tput rmul)

# Extract the name of this script.
cli=${0##*/}

# Define and create the setup log.
setuplog="setup.log"
> ${setuplog}

# Define the directory that will hold dependencies.
ecspert_deps=~/ecspert/deps
petsc_arch=arch-ecspert

# Define the target PETSc and SLEPc packages to download, if necessary.
petsc_package=petsc-3.20.2.tar.gz
slepc_package=slepc-3.20.1.tar.gz

# Set option defaults.
verbose=0
download_petsc=0
download_slepc=0

# This is the CLI's main help text.
show_help()
{
    echo "
${textbf}NAME${textnm}
        $cli - Set up the hybrid simulation and potential solver

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli${textnm} [${startul}OPTION${endul}]

${textbf}DESCRIPTION${textnm}
        This script will perform initial setup operations for ECSPERT.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print informative messages during the setup process.
        ${textbf}--download-petsc${textnm}
                Download and configure PETSc.
                This will install PETSc in ${ecspert_deps}, under arch-ecspert.
        ${textbf}--download-slepc${textnm}
                Download and configure SLEPc.
                This will install SLEPc in ${ecspert_deps}, under arch-ecspert.
        ${textbf}--dry-run${textnm}
                Display the sequence of commands to be executed.
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
    -o 'hv' \
    -l 'help,verbose' \
    -l 'download-petsc' \
    -l 'download-slepc' \
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
        '-v'|'--verbose')
            verbose=1
            shift
            continue
        ;;
        '--download-petsc')
            download_petsc=1
            shift
            continue
        ;;
        '--download-slepc')
            download_slepc=1
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
    echo "Did you misspell something?"
    echo
    echo "This script does not support passing arbitrary arguments for installing dependencies."
    echo "You may install individual dependencies on their own, then use the --with-<...> options"
    echo "to point to those installations."
    echo
    exit 1
fi

stage=

cleanup() {
    if [ "${stage}" != "${success}" ]; then
        if [ -n "${stage}" ]; then
            exitmsg="${stage} stage failed. See ${setuplog} for details."
        else
            exitmsg="Unknown error."
        fi
        echo 
        echo -e ${u_error} $exitmsg
        exit 1
    else
        echo 
        echo -e "${u_success} Success! ${u_success}"
        echo 
        echo "You may now try building the simulation or potential solver by running"
        echo "$ ./build.sh -p {pic,solver} [OPTIONS]"
    fi
}

trap cleanup EXIT

mark_stage() {
    stage="${1}"
    if [ $verbose == 1 ]; then
        echo "[${cli}] Executing ${stage} stage"
    fi
}

stage="installation"

# Adapted from https://stackoverflow.com/a/3352015. This will trim leading and
# trailing whitespace from a given string. For example: `trim "  a b  "` yields
# `"a b"`.
trim() {
    local var="$*"
    # leading whitespace
    var="${var#"${var%%[![:space:]]*}"}"
    # trailing whitespace
    var="${var%"${var##*[![:space:]]}"}"
    printf '%s' "$var"
}

# The function that will download and build external dependencies.
install_ext_deps() {
    if [ -z "${1}" ]; then
        echo "${u_error} Missing path to external-dependencies directory."
        exit 1
    fi
    local deps_dir="${1}"
    local tmp_name=tmp

    # Create and enter the top-level directory for external dependencies.
    mkdir -p ${deps_dir}
    pushd ${deps_dir} &> /dev/null

    # Create and enter the local subdirectory where this script will download
    # and build each package. Isolating the source code allows us to remove it
    # while leaving the build dependencies.
    mkdir ${tmp_name}
    pushd ${tmp_name} &> /dev/null

    # --> PETSc
    wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc_package
    tar -xvzf petsc_package
    pushd ${ecspert_deps}/petsc &> /dev/null
    export PETSC_ARCH=${petsc_arch} \
    ./configure \
    --download-cmake \
    --download-fblaslapack \
    --download-hdf5 \
    --download-scalapack \
    --download-mumps \
    --download-hypre \
    && make all check
    popd &> /dev/null
    export PETSC_DIR=${ecspert_deps}/petsc

    # --> SLEPc
    wget https://slepc.upv.es/download/distrib/slepc_package
    tar -xvzf slepc_package
    pushd ${ecspert_deps}/slepc &> /dev/null
    ./configure
    popd &> /dev/null
    export SLEPC_DIR=${ecspert_deps}/slepc

    # Exit and remove the temporary source-code directory.
    popd &> /dev/null
    /bin/rm -rf ${tmp_name}

    # Exit the top-level external-dependencies directory.
    popd &> /dev/null
}

# Design Questions
# - Should we reduce the options to a single --download-deps option?
# - Should we split install_ext_deps into a function for each dependency or
#   should we pass a flag to indicate which package(s) to install?
if [ -n "${download_petsc}" ]; then
    # TODO
fi

