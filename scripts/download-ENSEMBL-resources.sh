#!/usr/bin/env bash

###############################################################################
#
#   Download and extract genome resources from
#   ENSEMBL server (annotation in GTF and sequences in FASTA).
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 21-03-2020
#   LICENSE: Apache_2.0
#   USAGE:
#   bash download-ENSEMBL-resources.sh -s {hsa|mmu} -o {}
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
if [ -d $"${OUTDIR}" ]; then
    echo "Output directory already exists!"
    exit 1
fi

# create output directory
mkdir $"${OUTDIR}"

# download and extract FASTA and GTF files
case $SPECIES in
    hsa) # Homo sapiens
        URL="ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/"
        URL+="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        curl $URL \
            --output $"${OUTDIR}"/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        URL="ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/"
        URL+="Homo_sapiens.GRCh38.102.gtf.gz"
        curl $URL \
            --output $"${OUTDIR}"/Homo_sapiens.GRCh38.102.gtf.gz
        gunzip -c $"${OUTDIR}"/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
            1> $"${OUTDIR}"/Homo_sapiens.GRCh38.dna.primary_assembly.fa
        gunzip -c $"${OUTDIR}"/Homo_sapiens.GRCh38.102.gtf.gz \
            1> $"${OUTDIR}"/Homo_sapiens.GRCh38.102.gtf
        ;;
    mmu) # Mus musculus
        URL="ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/"
        URL+="Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
        curl $URL \
            --output $"${OUTDIR}"/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
        URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/"
        URL+="Mus_musculus.GRCm38.102.gtf.gz"
        curl $URL \
            --output $"${OUTDIR}"/Mus_musculus.GRCm38.102.gtf.gz
        gunzip -c $"${OUTDIR}"/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
            1> $"${OUTDIR}"/Mus_musculus.GRCm38.dna.primary_assembly.fa
        gunzip -c $"${OUTDIR}"/Mus_musculus.GRCm38.102.gtf.gz \
            1> $"${OUTDIR}"/Mus_musculus.GRCm38.102.gtf
        ;;
esac
