#!/usr/bin/env bash

###############################################################################
#
#   Download PolyAsite atlas (and select canonical chromosomes +MT)
#   for a given organism (supported: hsa/mmu).
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 21-03-2020
#   LICENSE: Apache_2.0
#   USAGE:
#   bash download-polyA-atlas.sh -s {hsa|mmu} -o {}
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
        -s|--species)
            SPECIES="$2"
            shift # past argument
            shift # past value
            ;;
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
if [ -z "$SPECIES" ] || [ -z "$OUTDIR" ]; then
    echo "Invalid arguments. Please provide --species and --output-directory"
    exit 1
fi

# check species option value
if [ ! "${SPECIES}" == "hsa" ] && [ ! "${SPECIES}" == "mmu" ]; then
    echo "Invalid species option. Currently supported are: 'hsa', 'mmu'"
    exit 1
fi

# exit if output directory exists
if [ -d "${OUTDIR}" ]; then
    echo "Output directory already exists!"
    exit 1
fi

# create directories
mkdir "${OUTDIR}"

# download and extract the Atlas
case $SPECIES in
    hsa) # Homo sapiens
        curl https://www.polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz \
            --output "${OUTDIR}"/atlas.clusters.2.0.GRCh38.96.bed.gz
        gunzip -c "${OUTDIR}"/atlas.clusters.2.0.GRCh38.96.bed.gz > "${OUTDIR}"/atlas.bed
        ;;
    mmu) # Mus musculus
        curl https://www.polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz \
            --output "${OUTDIR}"/atlas.clusters.2.0.GRCm38.96.bed.gz
        gunzip -c "${OUTDIR}"/atlas.clusters.2.0.GRCm38.96.bed.gz > "${OUTDIR}"/atlas.bed
        ;;
esac

# Filter canonical chromosomes (+MT)
# At this point we do not know what ENSEMBL annotation will be used in the pipeline
case $SPECIES in
    hsa) # Homo sapiens
        awk '$1 {if ($1 == "1") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "2") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "3") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "4") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "5") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "6") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "7") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "8") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "9") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "10") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "11") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "12") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "13") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "14") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "15") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "16") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "17") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "18") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "19") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "20") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "21") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "22") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "X") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "Y") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "MT") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        ;;
    mmu) # Mus musculus
        awk '$1 {if ($1 == "1") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "2") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "3") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "4") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "5") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "6") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "7") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "8") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "9") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "10") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "11") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "12") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "13") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "14") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "15") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "16") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "17") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "18") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "19") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "X") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "Y") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        awk '$1 {if ($1 == "MT") print}' "${OUTDIR}"/atlas.bed >> "${OUTDIR}"/filtered-atlas.bed
        ;;
esac

# Reformat original pasAtlas 2.0 to the standard MAPP atlas expected format
cat "${OUTDIR}"/filtered-atlas.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' > "${OUTDIR}"/clean-atlas.bed
