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
u_failure="\U1F62d"

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

# Determine this script's parent directory.
script_dir="$(dirname "$(readlink -f "${0}")")"

# Define the directory that will hold dependencies.
ecsprt_deps=${script_dir}/deps

# Define and create the setup log.
setuplog=${script_dir}/"setup.log"
> ${setuplog}

# Define the value to use for PETSC_ARCH with --download-petsc.
petsc_arch=arch-ecsprt

# Define the target PETSc and SLEPc packages to download, if necessary.
petsc_version="3.20.2"
slepc_version="3.20.1"
petsc_package=petsc-3.20.2.tar.gz
slepc_package=slepc-3.20.1.tar.gz

# Set option defaults.
verbose=0
download_petsc=0
download_slepc=0
petsc_dir=
slepc_dir=
petsc_arch=
create_dev=1

# This is the CLI's main help text.
show_help()
{
    echo "
${textbf}NAME${textnm}
        $cli - Set up the hybrid simulation and potential solver

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli${textnm} [${startul}OPTION${endul}]

${textbf}DESCRIPTION${textnm}
        This script will perform initial setup operations for ECSPRT.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print informative messages during the setup process.
        ${textbf}--download-petsc${textnm}
                Download and configure PETSc.
                This will install PETSc in ${ecsprt_deps}, under arch-ecsprt.
        ${textbf}--download-slepc${textnm}
                Download and configure SLEPc.
                This will install SLEPc in ${ecsprt_deps}, under arch-ecsprt.
        ${textbf}--with-petsc-dir DIR${textnm}
                Use DIR as PETSC_DIR when building PETSc and SLEPc.
                Ignored with --download-petsc.
        ${textbf}--with-slepc-dir DIR${textnm}
                Use DIR as SLEPC_DIR when building SLEPc.
                Ignored with --download-petsc.
        ${textbf}--with-petsc-arch ARCH${textnm}
                Use ARCH as PETSC_ARCH when building PETSc and SLEPc.
                Ignored with --download-petsc.
        ${textbf}--no-dev${textnm}
                Do not create the local dev/ subdirectory.
                By default, this script will create a directory called 
                ${script_dir}/dev.
                Git will automatically ignore its contents, making it useful for storing
                files that are specific to this copy of the repository (e.g., notes or 
                temporary scripts). This script will not overwrite dev/ if it exists.
"
}

# This will run if the CLI gets an unrecognized option.
report_bad_arg()
{
    printf "\nUnrecognized command: ${1}\n\n"
}

# Read command-line arguments.
TEMP=$(getopt \
    -n 'setup.sh' \
    -o 'hv' \
    -l 'help,verbose' \
    -l 'download-petsc,download-slepc' \
    -l 'with-petsc-dir:,with-slepc-dir:,with-petsc-arch:' \
    -l 'no-dev' \
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
        '--with-petsc-dir')
            export PETSC_DIR="${2}"
            shift 2
            continue
        ;;
        '--with-slepc-dir')
            export SLEPC_DIR="${2}"
            shift 2
            continue
        ;;
        '--with-petsc-arch')
            export PETSC_ARCH="${2}"
            shift 2
            continue
        ;;
        '--no-dev')
            create_dev=0
            shift 1
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
        echo -e ${u_failure} $exitmsg
        exit 1
    else
        echo 
        echo -e "${u_success} Success! ${u_success}"
        echo
        echo "You may want to set the following environment variables in your start-up file."
        echo "PETSC_DIR: \"${PETSC_DIR}\""
        echo "SLEPC_DIR: \"${SLEPC_DIR}\""
        echo 
        echo "For example, if you're running bash, add the following lines to ~/.bashrc"
        echo "export PETSC_DIR=\"${PETSC_DIR}\""
        echo "export SLEPC_DIR=\"${SLEPC_DIR}\""
        echo 
        echo "You should now try building the simulation or potential solver by running"
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

print_width() {
    local c=${1}

    printf "=%.0s" $(seq 1 $(tput cols))
    echo 
}

print_header() {
    print_width "="
    if [ -n "${@}" ]; then
        echo -e ${@}
        print_width "="
    fi
}

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

stage="installation"

build_petsc() {
    export PETSC_DIR=$(pwd)
    ./configure \
    --download-cmake \
    --download-fblaslapack \
    --download-hdf5 \
    --download-scalapack \
    --download-mumps \
    --download-hypre \
    &>> ${setuplog} \
    && make all check \
    &>> ${setuplog}
}

build_slepc() {
    export SLEPC_DIR=$(pwd)
    ./configure &>> ${setuplog} \
    && make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} \
    &>> ${setuplog} \
    && make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} check \
    &>> ${setuplog}
}

# The function that will download and build external dependencies.
install_ext_dep() {
    if [ -z "${1}" ]; then
        return 0
    fi

    local package_name="${1}"

    # Create and enter the top-level directory for external dependencies.
    mkdir -p ${ecsprt_deps}
    pushd ${ecsprt_deps} &> /dev/null

    # --> PETSc
    if [ ${package_name} == "petsc" ]; then
        print_header "Downloading PETSc"
        archive="petsc-${petsc_version}.tar.gz"
        wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${archive}
        print_header "Unpacking PETSc"
        tar -xvzf ${archive} &>> ${setuplog}
        ln -s "petsc-${petsc_version}" petsc
        pushd ${ecsprt_deps}/petsc &> /dev/null
        print_header "Building PETSc\nThis may take tens of minutes"
        build_petsc
        popd &> /dev/null
    fi

    # --> SLEPc
    if [ ${package_name} == "slepc" ]; then
        print_header "Downloading SLEPc"
        archive="slepc-${slepc_version}.tar.gz"
        wget https://slepc.upv.es/download/distrib/${archive}
        print_header "Unpacking SLEPc"
        tar -xvzf ${archive} &>> ${setuplog}
        ln -s "slepc-${slepc_version}" slepc
        pushd ${ecsprt_deps}/slepc &> /dev/null
        print_header "Building SLEPc\nThis may take a few minutes"
        build_slepc
        popd &> /dev/null
    fi

    # Exit the top-level external-dependencies directory.
    popd &> /dev/null
}

# Download and install PETSc, if requested.
if [ $download_petsc == 1 ]; then
    install_ext_dep "petsc"
fi

# Download and install SLEPc, if requested.
if [ $download_slepc == 1 ]; then
    install_ext_dep "slepc"
fi

# Offer to create git-ignored dev directory.
if [ $create_dev == 1 ]; then
    mkdir -p dev
fi

stage=$success

