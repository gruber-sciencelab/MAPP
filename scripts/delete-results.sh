#!/usr/bin/env bash

###############################################################################
#
#   Delete all output after a MAPP run.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 27-11-2020
#   LICENSE: Apache_2.0
#   USAGE: bash scripts/delete-results.sh
#
###############################################################################

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

for dir in "$SCRIPTDIR"/../modules/*/
do
    rm -rf "$dir"/output
done

rm -f "$SCRIPTDIR"/../summary.tar.gz
rm -rf "$SCRIPTDIR"/../summary "$SCRIPTDIR"/../logs
