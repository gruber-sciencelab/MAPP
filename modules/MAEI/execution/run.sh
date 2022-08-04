###############################################################################
#
#   Main execution script for the MAE workflow
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 29-11-2021
#   LICENSE: Apache_2.0
#
###############################################################################

# exit at a first command that exits with a !=0 status
set -eo pipefail

user_dir=$PWD
pipeline_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd "$pipeline_dir"

# in case of SIGINT go back to the user dir
goback () {
    rc=$?
    cd "$user_dir"
    echo "Exit status: $rc"
}
trap goback SIGINT

###############################################################################
# parse command line arguments
###############################################################################

# no args
if [ $# -eq 0 ]
then
    echo "Please execute $ run.sh --help to learn how to call this script."
fi

# --help msg
if printf '%s\n' "$@" | grep -q '^--help$'; then
    echo "This is the main script to call the MAE workflow."
    echo "Available options:"
    echo ""
    echo "  -c/--configfile {XXX}"
    echo "  Path to the snakemake config file."
    echo ""
    echo "  -e/--environment {local/slurm/sge} (REQUIRED)"
    echo "  Environment to execute the workflow in:"
    echo "  * local = execution on the local machine."
    echo "  * slurm = slurm cluster support."
    echo "  * sge = (sun) grid engine cluster support."
    echo ""
    echo "  -t/--technology {conda/singularity} (REQUIRED)"
    echo "  Technology for reproducible research:"
    echo "  * conda = use conda environments throughout the workflow"
    echo "  * singularity = use singularity containers throughout the workflow"
    echo ""
    echo "  -b/--bind {XXX,YYY} (OPTIONAL)"
    echo "  For workflow execution in a cluster env. with singularity tech:"
    echo "  additional ABSOLUTE paths that need to be accessible from the containers."
    echo ""
    echo "  -g/--graph {rulegraph/dag} {XXX} (OPTIONAL)"
    echo "  Do not call the execution."
    echo "  Instead - generate Snakemake graps of the workflow in XXX file (SVG)."
    echo "  * rulegraph = create a rule graph"
    echo "  * dag = create a directed acyclic graph"
    exit 0
fi

# parse args
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -c|--configfile)
            CONFIGFILE="$2"
            shift # past argument
            shift # past value
            ;;
        -e|--environment)
            ENV="$2"
            shift # past argument
            shift # past value
            ;;
        -t|--technology)
            TECH="$2"
            shift # past argument
            shift # past value
            ;;
        -b|--bind)
            BINDPATHS="$2"
            shift # past argument
            shift # past value
            ;;
        -g|--graph)
            GRAPH="$2"
            shift # past argument
            GRAPHPATH="$2"
            shift # past value
            shift # past value
            ;;
        *) # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

###############################################################################
# MAIN
###############################################################################

# configfile is always required
if [ -z "$CONFIGFILE" ]; then
    echo "Invalid arguments. Please provide --configfile"
    exit 1
fi
if ! [[ "$CONFIGFILE" = /* ]]; then
    CONFIGFILE=$user_dir/$CONFIGFILE
fi
if ! [ -f "$CONFIGFILE" ]; then
    echo "Invalid arguments. Configuration file does not exist"
    exit 1
fi

# if this is just a graph run we don't need other parameters
if [ -n "$GRAPH" ]; then

    if [ "$GRAPH" != "rulegraph" ] && [ "$GRAPH" != "dag" ]; then
        echo "Invalid arguments. --graph = {rulegraph/dag}"
        exit 1
    fi

    if [ "$GRAPH" == "rulegraph" ]; then
        snakemake \
            --snakefile="../Snakefile" \
            --configfile="$CONFIGFILE" \
            --rulegraph \
            | dot -Tsvg \
            > "$user_dir/$GRAPHPATH"
    fi

    if [ "$GRAPH" == "dag" ]; then
        snakemake \
            --snakefile="../Snakefile" \
            --configfile="$CONFIGFILE" \
            --dag \
            | dot -Tsvg \
            > "$user_dir/$GRAPHPATH"
    fi

    exit 0
fi

# if this is not a graph run - check other parameters
if [ -z "$ENV" ]; then
    echo "Invalid arguments. Please provide --environment"
    exit 1
fi
if [ "$ENV" != "local" ] && [ "$ENV" != "slurm" ] && [ "$ENV" != "sge" ]; then
    echo "Invalid arguments. --environment = {local/slurm/sge}"
    exit 1
fi
if [ -z "$TECH" ]; then
    echo "Invalid arguments. Please provide --technology"
    exit 1
fi
if [ "$TECH" != "conda" ] && [ "$TECH" != "singularity" ]; then
    echo "Invalid arguments. --technology = {conda/singularity}"
    exit 1
fi

# select a proper smk profile based on the command line args
case "$ENV$TECH" in
    localconda)
        snakemake \
            --configfile="$CONFIGFILE" \
            --profile="../profiles/local-conda"
        ;;
    localsingularity)
        snakemake \
            --configfile="$CONFIGFILE" \
            --profile="../profiles/local-singularity" \
            --singularity-args "--no-home --bind ${PWD}/..,$BINDPATHS"
        ;;
    slurmconda)
        snakemake \
            --configfile="$CONFIGFILE" \
            --profile="../profiles/slurm-conda"
        ;;
    slurmsingularity)
        snakemake \
            --configfile="$CONFIGFILE" \
            --profile="../profiles/slurm-singularity" \
            --singularity-args "--no-home --bind ${PWD}/..,$BINDPATHS"
        ;;
    sgeconda)
        snakemake \
            --configfile="$CONFIGFILE" \
            --profile="../profiles/sge-conda"
        ;;
    sgesingularity)
        echo "Singularity technology is not supported for Grid Engine yet."
        ;;
esac
