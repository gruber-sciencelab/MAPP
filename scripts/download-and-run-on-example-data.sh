#!/usr/bin/env bash

###############################################################################
#
#   Example MAPP workflow run:
#   * download test data
#   * set up the environment
#   * prepare configuration files
#   * initiate local run in conda envs
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 13-10-2021
#   LICENSE: Apache_2.0
#
###############################################################################

# parse the optional number of cores
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -n|--cores)
            CORES="$2"
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

if [ -z "$CORES" ]; then
    CORES=1
fi

echo "### Downloading test dataset ###"
wget https://zenodo.org/record/6941629/files/MAPP_test_data.tar.gz

echo "### Extracting test dataset ###"
tar -xzf MAPP_test_data.tar.gz
rm -rf MAPP_test_data.tar.gz

echo "### Adjusting configuration template ###"
MAPP_directory=$(pwd)

# Linux version
if [ "$(uname)" == "Linux" ]; then

    # Specify absolute paths to all the test data:
    sed -i '23s|""|"'"$MAPP_directory"'"|'  MAPP_test_data/config_template.yml
    genomic_annotation=${MAPP_directory}/MAPP_test_data/Homo_sapiens.GRCh38.87.chr21.gtf
    sed -i '26s|""|"'"$genomic_annotation"'"|'  MAPP_test_data/config_template.yml
    genomic_sequence=${MAPP_directory}/MAPP_test_data/Homo_sapiens.GRCh38.dna.chromosome.21.fa
    sed -i '35s|""|"'"$genomic_sequence"'"|'  MAPP_test_data/config_template.yml
    genomic_index=${MAPP_directory}/MAPP_test_data/GRCh38_STAR_index_hsa105
    sed -i '41s|""|"'"$genomic_index"'"|'  MAPP_test_data/config_template.yml
    analysis_design_table=${MAPP_directory}/MAPP_test_data/design_table.tsv
    sed -i '44s|""|"'"$analysis_design_table"'"|'  MAPP_test_data/config_template.yml
    PWM_directory=${MAPP_directory}/MAPP_test_data/PWMS
    sed -i '48s|""|"'"$PWM_directory"'"|'  MAPP_test_data/config_template.yml
    seqlogo_directory=${MAPP_directory}/MAPP_test_data/SEQLOGOS
    sed -i '52s|""|"'"$seqlogo_directory"'"|'  MAPP_test_data/config_template.yml
    PAS_atlas=${MAPP_directory}/MAPP_test_data/polyAsiteAtlas2.chr21.bed
    sed -i '55s|""|"'"$PAS_atlas"'"|'  MAPP_test_data/config_template.yml

    #echo "### Activating conda env for MAPP workflow ###"
    source "$CONDA_PREFIX"/etc/profile.d/conda.sh
    conda activate mapp

fi

# Mac version
if [ "$(uname)" == "Darwin" ]; then

    # Specify absolute paths to all the test data:
    sed -i '' -e '23s|""|"'"$MAPP_directory"'"|'  MAPP_test_data/config_template.yml
    genomic_annotation=${MAPP_directory}/MAPP_test_data/Homo_sapiens.GRCh38.87.chr21.gtf
    sed -i '' -e '26s|""|"'"$genomic_annotation"'"|'  MAPP_test_data/config_template.yml
    genomic_sequence=${MAPP_directory}/MAPP_test_data/Homo_sapiens.GRCh38.dna.chromosome.21.fa
    sed -i '' -e '35s|""|"'"$genomic_sequence"'"|'  MAPP_test_data/config_template.yml
    genomic_index=${MAPP_directory}/MAPP_test_data/GRCh38_STAR_index_hsa105
    sed -i '' -e '41s|""|"'"$genomic_index"'"|'  MAPP_test_data/config_template.yml
    analysis_design_table=${MAPP_directory}/MAPP_test_data/design_table.tsv
    sed -i '' -e '44s|""|"'"$analysis_design_table"'"|'  MAPP_test_data/config_template.yml
    PWM_directory=${MAPP_directory}/MAPP_test_data/PWMS
    sed -i '' -e '48s|""|"'"$PWM_directory"'"|'  MAPP_test_data/config_template.yml
    seqlogo_directory=${MAPP_directory}/MAPP_test_data/SEQLOGOS
    sed -i '' -e '52s|""|"'"$seqlogo_directory"'"|'  MAPP_test_data/config_template.yml
    PAS_atlas=${MAPP_directory}/MAPP_test_data/polyAsiteAtlas2.chr21.bed
    sed -i '' -e '55s|""|"'"$PAS_atlas"'"|'  MAPP_test_data/config_template.yml

    #echo "### Activating conda env for MAPP workflow ###"
    source "$CONDA_PREFIX"/etc/profile.d/conda.sh
    conda activate mapp

fi

echo "### Setting up configuration file for MAPP workflow ###"
python scripts/create-main-config-file.py \
    --config-template MAPP_test_data/config_template.yml \
    --pipeline-configfile configs/config.yml

#echo "### Initiating MAPP workflow ###"
bash execution/run.sh \
    --configfile configs/config.yml \
    --environment local  \
    --technology conda \
    --cores "$CORES"
