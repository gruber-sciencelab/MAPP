#!/usr/bin/env bash

###############################################################################
#
#   Create snakemake report for the whole MAPP run.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 27-11-2020
#   LICENSE: Apache_2.0
#   USAGE:
#   bash create-snakemake-report.sh -r {}
#
###############################################################################

###############################################################################
# parse command line arguments
###############################################################################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -r|--report-path)
            REPPATH="$2"
            shift # past argument
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

# check arguments
if [ -z "$REPPATH" ]; then
    echo "Invalid arguments. Please provide --report-path"
    exit 1
fi

cleanup () {
    rc=$?
    cd "$user_dir"
    echo "Exit status: $rc"
}
trap cleanup SIGINT

set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands

# set paths
user_dir=$PWD
execution_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"/../execution
if [[ "$REPPATH" != /* ]]; then
    REPPATH=$user_dir/"$REPPATH"
fi

# generate the snakemake report
cd "$execution_dir"
snakemake \
    --snakefile="../Snakefile" \
    --configfile="../configs/config.yml" \
    --report=$"${REPPATH}"
