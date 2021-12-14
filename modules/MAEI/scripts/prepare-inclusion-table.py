"""
##############################################################################
#
#   Filter out values in the exon inclusion table for the subsequent
#   MARA model. Also, change the exon coordinates system from GTF to BED.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 09-04-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
from itertools import compress
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
        "--inclusion-table",
        dest="inclusion_table",
        required=True,
        help="Path to the exon inclusion table.",
    )
    parser.add_argument(
        "--filtered-inclusion-table",
        dest="filtered_inclusion_table",
        required=True,
        help="Path for the filtered exon inclusion table.",
    )
    parser.add_argument(
        "--minimal-expression",
        dest="min_exp",
        required=True,
        help="Minimal expression of all transcripts supporting an AS exon (for all samples).",
    )
    return parser


##############################################################################


def convert_format(coord_str):
    """
    Subtracts one from the beg coordinate.
    """
    return (
        coord_str.split(":")[0]
        + ":"
        + str(int(coord_str.split(":")[1].split("-")[0]) - 1)
        + "-"
        + coord_str.split("-")[1].split(":")[0]
        + ":"
        + coord_str[-1]
    )


def main():
    """Main body of the script."""

    # generate a set of exons from inclusion table
    inclusion_table_df = pd.read_csv(options.inclusion_table, sep="\t")
    # the coordinates will later have to match those in bed
    # BED files are 0 based, half open at the end. GTFs are 1-based, ends closed
    inclusion_table_df["coordinates"] = inclusion_table_df.apply(
        lambda row: convert_format(row.coordinates), axis=1
    )
    inclusion_table_exons = list(inclusion_table_df["coordinates"])
    keep_exon_marker = [True] * len(inclusion_table_exons)

    # separate the df into included/total tables
    included_cols = [
        c for c in inclusion_table_df.columns.values if c.find("_included") > -1
    ]
    total_cols = [c for c in inclusion_table_df.columns.values if c.find("_total") > -1]
    included_df = inclusion_table_df[included_cols].copy()
    included_df.columns = list(range(len(included_df.columns)))
    total_df = inclusion_table_df[total_cols].copy()
    total_df.columns = list(range(len(total_df.columns)))

    # filter 1: apply minimal expression on all transcripts supporting given AS exon (for all samples)
    min_exp = float(options.min_exp)
    below_cutoff_marker_list = list((total_df < min_exp).any(axis=1))
    for i, x in enumerate(below_cutoff_marker_list):
        if x:
            keep_exon_marker[i] = False

    # filter 2: Filter exons which are not included in any sample (not informative for the model)
    inclusion_epsilon = 0.00001
    below_epsilon_marker_list = list((included_df < inclusion_epsilon).all(axis=1))
    for i, x in enumerate(below_epsilon_marker_list):
        if x:
            keep_exon_marker[i] = False

    # filter 3: Filter exons which are fully included in all samples (not informative for the model)
    delta_epsilon = 0.00001
    diff_df = total_df - included_df
    below_delta_epsilon_marker_list = list((diff_df < delta_epsilon).all(axis=1))
    for i, x in enumerate(below_delta_epsilon_marker_list):
        if x:
            keep_exon_marker[i] = False

    # apply the filtering
    inclusion_table_exons = list(compress(inclusion_table_exons, keep_exon_marker))

    # filter the inclusion table
    inclusion_table_df = inclusion_table_df[
        inclusion_table_df["coordinates"].isin(inclusion_table_exons)
    ]
    inclusion_table_df = inclusion_table_df.set_index("coordinates").sort_index()
    inclusion_table_df.to_csv(options.filtered_inclusion_table, sep="\t")


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
