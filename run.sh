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

# Declare run-related directories and files.
cwd=$(pwd)
srcdir=${cwd}/src
bindir=${cwd}/bin
logname="run.log"
runlog=

# Set option defaults.
np=1
verbose=0
debug=0
rundir='.'
outdir=
options=
optlist=
reqopts=0
args=

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
        $cli - Run the hybrid simulation or potential solver

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli${textnm} -p {pic,solver} [${startul}OPTIONS${endul}] [args ${startul}ARGS${endul}]

${textbf}DESCRIPTION${textnm}
        This script will run the target program. All arguments after 'args' will pass through.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-p${textnm}, ${textbf}--program${textnm} ${startul}PROG${endul}
                The target program (required).
        ${textbf}--debug${textnm}
                Run ${startul}PROG${endul} in debugging mode.
        ${textbf}-n${textnm}, ${textbf}--nproc${textnm} ${startul}N${endul}
                Number of processors to use.
                (default: ${np})
        ${textbf}-r${textnm}, ${textbf}--rundir${textnm} ${startul}RUNDIR${endul}
                Name of the parent directory that should contain the 
                run-specific directory. See ${startul}OUTDIR${endul}.
                (default: ${rundir})
        ${textbf}-o${textnm}, ${textbf}--outdir${textnm} ${startul}OUTDIR${endul}
                Name of the run-specific directory, within ${startul}RUNDIR${endul}, that will
                contain the output of this run.
                The default name is generated from the current date and time.
        ${textbf}--options${textnm} ${startul}FILE${endul}
                Full name of the file that contains runtime options.
                You may pass this argument more than once to chain options files.
                (default: none)
        ${textbf}--options-list${textnm} ${startul}FILE${endul}
                Full name of a file that contains names of options files to merge into a 
                single file of runtime options.
                (default: none)
        ${textbf}--require-options${textnm}
                If true, raise an error if an options file does not exist or is not a regular file.
                (default: false)
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print informative messages during the execution process.
"
}

# This will run if the CLI gets an unrecognized option.
report_bad_arg()
{
    printf "\nUnrecognized command: ${1}\n\n"
}

usercmnd="$0 $@"

# Read command-line arguments.
TEMP=$(getopt \
    -n 'run.sh' \
    -o '+hvp:n:r:o:' \
    -l 'help,verbose' \
    -l 'program:' \
    -l 'debug' \
    -l 'nprocs:' \
    -l 'rundir:' \
    -l 'outdir:' \
    -l 'options:' \
    -l 'options-list:' \
    -l 'require-options' \
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
        '--debug')
            debug=1
            shift
            continue
        ;;
        '-n'|'--nproc')
            np="${2}"
            shift 2
            continue
        ;;
        '-r'|'--rundir')
            rundir="${2}"
            shift 2
            continue
        ;;
        '-o'|'--outdir')
            outdir="${2}"
            shift 2
            continue
        ;;
        '--options')
            current=$(realpath -m "${2}")
            if [ -z "${options}" ]; then
                options="${current}"
            else
                options="${options} ${current}"
            fi
            shift 2
            continue
        ;;
        '--options-list')
            filepath=$(realpath -m "${2}")
            parts=($(cat "${filepath}"))
            optlist=($(for part in ${parts[@]}; do echo $(realpath -m "${part}"); done))
            shift 2
            continue
        ;;
        '--require-options')
            reqopts=1
            shift
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

# Convert extra arguments to program arguments.
extra="$@"
args=${extra//args}

stage=

cleanup() {
    if [ "${stage}" != "${success}" ]; then
        if [ $verbose == 1 ]; then
            if [ -n "${stage}" ]; then
                exitmsg="${stage} stage failed. See ${runlog} for details."
            else
                exitmsg="Unknown error."
            fi
            echo 
            echo -e ${u_failure} $exitmsg
        fi
        exit 1
    else
        if [ $verbose == 1 ]; then
            echo 
            echo -e "${u_success} Success! ${u_success}"
            echo 
            echo "Output is in $(realpath -L ${dstdir})"
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

# Define a temporary run log path, in case setup fails.
runlog=${rundir}/${logname}

# Mark this stage.
mark_stage "setup"

# Set the local name of the output directory.
if [ -z "${outdir}" ]; then
    dstdir=${rundir}/$(date +'%Y-%m-%d-%H%M%S')
else
    dstdir=${rundir}/"${outdir}"
fi

# Create the full output directory.
mkdir -p ${dstdir}

# Redefine the run log path.
runlog=$(realpath -m ${dstdir}/${logname})
> ${runlog}

# Refuse to run without a target program.
if [ -z "${prog}" ]; then
    if [ $verbose == 1 ]; then
        show_help
        echo
        echo "ERROR: Missing target program." &> ${runlog}
    fi
    exit 1
fi

# Echo the executed command.
if [ $verbose == 1 ]; then
    echo "[${cli}] command: ${usercmnd}" &> ${runlog}
fi

# Mark this stage.
mark_stage "symlink"

# Create a symlink to this run in the directory of runs.
pushd ${rundir} &> /dev/null
rm -f latest
ln -s ${dstdir} latest
popd &> /dev/null

# Move to the output directory.
pushd ${dstdir} &> /dev/null

# Mark this stage.
mark_stage "options"

# Create the runtime options database.
tmpopts=.options.txt
touch "${tmpopts}"
if [ -n "${optlist}" ]; then
    for path in ${optlist[@]}; do
        if [ -f "${path}" ]; then
            if [ $verbose == 1 ]; then
                echo "[${cli}] Adding options from ${path}" &>> ${runlog}
            fi
            cat "${path}" >> "${tmpopts}"
            echo >> "${tmpopts}"
        else
            message="Cannot add options from ${path}: file does not exist or is not a regular file"
            if [ $reqopts == 1 ]; then
                if [ $verbose == 1 ]; then
                    echo "[${cli}] ERROR: ${message}" &>> ${runlog}
                fi
                exit 1
            fi
            if [ $verbose == 1 ]; then
                echo "[${cli}] WARNING: ${message}" &>> ${runlog}
            fi
        fi
    done
fi
if [ -n "${options}" ]; then
    for current in ${options[@]}; do
        if [ -f "${current}" ]; then
            if [ $verbose == 1 ]; then
                echo "[${cli}] Adding options from ${current}" &>> ${runlog}
            fi
            (/bin/cat "${current}" >> "${tmpopts}") &>> ${runlog}
        else
            message="Cannot add options from ${current}: file does not exist or is not a regular file"
            if [ $reqopts == 1 ]; then
                if [ $verbose == 1 ]; then
                    echo "[${cli}] ERROR: ${message}" &>> ${runlog}
                fi
                exit 1
            fi
            if [ $verbose == 1 ]; then
                echo "[${cli}] WARNING: ${message}" &>> ${runlog}
            fi
        fi
    done
fi
if [ -s "${tmpopts}" ]; then
    /bin/cp "${tmpopts}" petsc.ini &>> ${runlog}
fi
rm ${tmpopts}

# Mark this stage.
mark_stage "run"

# Run the program.
if [ ${debug} == 1 ]; then
    if [ ${np} -gt 1 ]; then
        if [ $verbose == 1 ]; then
            echo -e "ERROR: Multi-processor debugging is currently disabled." &>> ${runlog}
        fi
        exit 1
    else
        gdb --args ${bindir}/${prog} $extra
    fi
else
    mpiexec -n ${np} ${bindir}/${prog} $extra &>> ${runlog}
fi

# Signal success to the clean-up function.
stage=${success}

popd &> /dev/null
