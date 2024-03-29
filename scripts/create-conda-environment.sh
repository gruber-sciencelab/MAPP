#!/usr/bin/env bash

###############################################################################
#
#   Create conda virtual environment for MAPP
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 20-03-2020
#   LICENSE: Apache_2.0
#   USAGE: bash create-conda-environment.sh
#
###############################################################################

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
mamba env create --file "$CWD"/../env/environment.yml
