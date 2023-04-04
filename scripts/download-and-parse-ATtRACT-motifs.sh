#!/usr/bin/env bash

###############################################################################
#
#   Download ATtRACT database and extract a clean set of hsa PWMs
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 21-03-2020
#   LICENSE: Apache_2.0
#   USAGE:
#   bash download-and-parse-ATtRACT-motifs.sh -o {}
#
###############################################################################

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

###############################################################################
# parse command line arguments
###############################################################################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -o|--output-directory)
            OUTDIR="$2"
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
if [ -z "$OUTDIR" ]; then
    echo "Invalid arguments. Please provide --output-directory"
    exit 1
fi

# exit if output directory exists
if [ -d "${OUTDIR}" ]; then
    echo "Output directory already exists!"
    exit 1
fi

# create output directory
mkdir "${OUTDIR}"

# download and unzip ATtRACT database
curl https://attract.cnic.es/attract/static/ATtRACT.zip \
    --output "${OUTDIR}"/ATtRACT.zip \
    2> /dev/null
unzip "${OUTDIR}"/ATtRACT.zip -d "${OUTDIR}" 1> /dev/null

# clean the set of motifs (custom filtering)
python scripts/clean-ATtRACT-motifs.py \
--dbpath "${OUTDIR}"/ATtRACT_db.txt \
--pwmpath "${OUTDIR}"/pwm.txt \
--outdir "${OUTDIR}"/motifs

# create sequence logos for every motif
mkdir "${OUTDIR}"/seqlogos
Rscript "$CWD"/create-seqlogos.r \
    --PWM_dir "${OUTDIR}"/motifs \
    --output_dir "${OUTDIR}"/seqlogos
