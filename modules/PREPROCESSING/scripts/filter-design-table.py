"""
##############################################################################
#
#   Filter design table based on the quality of RNA-Seq samples.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 02-05-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import os
from datetime import datetime
import numpy as np
import pandas as pd


def parse_arguments():
    """Parser of the command-line arguments."""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-v",
        "--verbosity",
        dest="verbosity",
        choices=("DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"),
        default="ERROR",
        help="Verbosity/Log level. Defaults to ERROR",
    )
    parser.add_argument(
        "-l", "--logfile", dest="logfile", help="Store log to this file."
    )
    parser.add_argument(
        "--mapping-quality-scores",
        dest="mapping_quality_table",
        required=True,
        help="Merged RNA-SeQC results per sample in a TSV format.",
    )
    parser.add_argument(
        "--TIN-scores",
        dest="TIN",
        required=True,
        help="Median TIN scores per sample in a TSV format.",
    )
    parser.add_argument(
        "--design-table-in",
        dest="input",
        required=True,
        help="Original design table in a TSV format.",
    )
    parser.add_argument(
        "--design-table-out",
        dest="output",
        required=True,
        help="Filtered design table in a TSV format.",
    )
    parser.add_argument(
        "--min-median-TIN-cutoff",
        dest="TIN_cutoff",
        required=True,
        help="Minimal value for the median TIN score per sample.",
    )
    parser.add_argument(
        "--min-mapping-rate",
        dest="rnaseqc_min_mapping_rate",
        required=True,
        help="RNA-SeQC: minimal value for the Mapping Rate.",
    )
    parser.add_argument(
        "--min-unique-rate-of-mapped",
        dest="rnaseqc_min_unique_rate_of_mapped",
        required=True,
        help="RNA-SeQC: minimal value for the Unique Rate of Mapped.",
    )
    parser.add_argument(
        "--min-high-quality-rate",
        dest="rnaseqc_min_high_quality_rate",
        required=True,
        help="RNA-SeQC: minimal value for the High Quality Rate.",
    )
    parser.add_argument(
        "--max-intergenic-rate",
        dest="rnaseqc_max_intergenic_rate",
        required=True,
        help="RNA-SeQC: maximal value for the Intergenic Rate.",
    )
    parser.add_argument(
        "--max-rRNA-rate",
        dest="rnaseqc_max_rRNA_rate",
        required=True,
        help="RNA-SeQC: maximal value for the rRNA Rate.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    samples_to_remove = []

    # read the original design table
    design_table_in = pd.read_csv(options.input, sep="\t", index_col=0)

    # quality filtering based on TIN score calculation
    TIN = pd.read_csv(options.TIN, sep="\t", index_col=0, header=None)
    TIN.index.rename("sample", inplace=True)
    TIN.columns = ["median_TIN_score"]
    for sample_ID in TIN.index.values:
        if float(TIN.at[sample_ID, "median_TIN_score"]) < float(options.TIN_cutoff):
            samples_to_remove.append(sample_ID)

    # quality filtering based on the mapping scores (RNA-SeQC)
    rnaseqc_table = pd.read_csv(options.mapping_quality_table, sep="\t", index_col=0)

    # RNA-SeQC quality filtering
    for sample_ID in rnaseqc_table.columns.values:
        if (
            rnaseqc_table[sample_ID]["Mapping Rate"]
            < float(options.rnaseqc_min_mapping_rate)
            or rnaseqc_table[sample_ID]["Unique Rate of Mapped"]
            < float(options.rnaseqc_min_unique_rate_of_mapped)
            or rnaseqc_table[sample_ID]["High Quality Rate"]
            < float(options.rnaseqc_min_high_quality_rate)
            or rnaseqc_table[sample_ID]["Intergenic Rate"]
            > float(options.rnaseqc_max_intergenic_rate)
            or rnaseqc_table[sample_ID]["rRNA Rate"]
            > float(options.rnaseqc_max_rRNA_rate)
        ):
            samples_to_remove.append(sample_ID)

    # filter the design table
    now = datetime.now()
    samples_to_remove = list(set(samples_to_remove))
    samples_to_keep = [
        s for s in design_table_in.index.values if s not in samples_to_remove
    ]
    design_table_out = design_table_in.loc[samples_to_keep].copy()
    design_table_out.to_csv(options.output, sep="\t")

    # print warnings to standard error stream
    datetime_string = now.strftime("%d/%b/%Y %H:%M:%S")
    for s in samples_to_remove:
        sys.stderr.write("[" + datetime_string + "] WARNING - ")
        sys.stderr.write("Sample: " + s + " has been removed due to low quality." + os.linesep)

    # the filtered table is already saved, now check if there are any samples left after qc
    # if not: raise exception to shutdown the whole pipeline (keeping the output file)
    if len(samples_to_keep) < 2:
        errmsg = "Less than two samples passed the quality control step. Further processing is not possible."
        raise RuntimeWarning(errmsg)


##############################################################################

if __name__ == "__main__":

    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        formatter = logging.Formatter(
            fmt="[%(asctime)s] %(levelname)s - %(message)s",
            datefmt="%d-%b-%Y %H:%M:%S",
        )
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger = logging.getLogger("logger")
        logger.setLevel(logging.getLevelName(options.verbosity))
        logger.addHandler(console_handler)
        if options.logfile is not None:
            logfile_handler = logging.handlers.RotatingFileHandler(
                options.logfile, maxBytes=50000, backupCount=2
            )
            logfile_handler.setFormatter(formatter)
            logger.addHandler(logfile_handler)

        # execute the body of the script
        start_time = time.time()
        logger.info("Starting script")
        main()
        seconds = time.time() - start_time

        # log the execution time
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        logger.info(
            "Successfully finished in {hours}h:{minutes}m:{seconds}s",
            hours=int(hours),
            minutes=int(minutes),
            seconds=int(seconds) if seconds > 1.0 else 1,
        )
    # log the exception in case it happens
    except Exception as e:
        logger.exception(str(e))
        raise e
