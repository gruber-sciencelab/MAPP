#!/usr/bin/env bash

###############################################################################
#
#   Clean up large files and directories after MAPP run.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 27-11-2020
#   LICENSE: Apache_2.0
#   USAGE: bash scripts/cleanup-results.sh
#
###############################################################################

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

# PREPROCESSING:
# * reads after adapter removal
# * reads after tail removal
# * STAR-related data (including alignments)
rm -rf "$CWD"/../modules/PREPROCESSING/output/adapter_trimmed
rm -rf "$CWD"/../modules/PREPROCESSING/output/tail_trimmed
for dir in "$CWD"/../modules/PREPROCESSING/output/alignments/*/
do
    prefix=$dir
    rm -rf "$prefix"/*._STARgenome
    rm -rf "$prefix"/*._STARpass1
    rm -f "$prefix"/*.Aligned.out.bam
    rm -f "$prefix"/*.Aligned.out.sorted.bam
    rm -f "$prefix"/*.Aligned.out.sorted.bam.bai
    rm -f "$prefix"/*.Aligned.toTranscriptome.out.bam
done

# CREATE_SITECOUNT_MATRICES:
# * sitecount matrices (and their links) from each window
for dir in "$CWD"/../modules/CREATE_SITECOUNT_MATRICES/output/3ss/*/
do
    prefix=$dir
    rm -f "$prefix"/matrix.tsv
    rm -f "$prefix"/kmers_matrix.tsv
    rm -f "$prefix"/pwms_matrix.tsv
done
for dir in "$CWD"/../modules/CREATE_SITECOUNT_MATRICES/output/5ss/*/
do
    prefix=$dir
    rm -f "$prefix"/matrix.tsv
    rm -f "$prefix"/kmers_matrix.tsv
    rm -f "$prefix"/pwms_matrix.tsv
done
for dir in "$CWD"/../modules/CREATE_SITECOUNT_MATRICES/output/pas/*/
do
    prefix=$dir
    rm -f "$prefix"/matrix.tsv
    rm -f "$prefix"/kmers_matrix.tsv
    rm -f "$prefix"/pwms_matrix.tsv
done

# PAQR:
# * remove pas coverages files
rm -rf "$CWD"/../modules/PAQR/output/pas_coverages

# QUANTIFICATION:
# * transcriptome sequences
rm -f "$CWD"/../modules/QUANTIFICATION/output/transcriptome.fasta
