"""
##############################################################################
#
#   Filter tandem pas according to the expression
#   (keep tandem pas expressed in all samples)
#   (+ sort the output tables with pas expression & position)
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 11-11-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
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
        "--normalized-expression",
        dest="normalized_expression",
        required=True,
        help="Path to the TSV file with normalized expression values.",
    )
    parser.add_argument(
        "--pas-positions",
        dest="pas_positions",
        required=True,
        help="Path to the file with pas positions within terminal exons.",
    )
    parser.add_argument(
        "--filtered-expression",
        dest="filtered_expression",
        required=True,
        help="Path for the output file with pas expression.",
    )
    parser.add_argument(
        "--filtered-positions",
        dest="filtered_positions",
        required=True,
        help="Path for the output file with pas positions.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    positions = pd.read_csv(options.pas_positions, sep="\t", index_col=0, header=None)
    expression = pd.read_csv(options.normalized_expression, sep="\t")
    samples = expression.columns.values[10:]

    # all sites in this table are significant in at least one sample
    # (according to the previous processing)
    expression = expression.replace(to_replace=-1.0, value=np.nan)
    filtered_exp_df = expression.dropna(how="any")

    # remove pas which are single per a terminal exon
    pas_blacklist = []
    for te, group in filtered_exp_df.groupby("exon"):
        if len(group) == 1:
            pas_blacklist.append(list(group["pas"])[0])
    filtered_exp_df = filtered_exp_df[~filtered_exp_df["pas"].isin(pas_blacklist)]

    # sort both expression and position tables accordingly and save
    filtered_exp_df.sort_values(by="pas", inplace=True)
    filtered_exp_df.to_csv(options.filtered_expression, sep="\t", index=False)
    pas_whitelist = list(filtered_exp_df["pas"])
    positions.index.name = "pas"
    positions = positions.loc[pas_whitelist]
    positions.sort_index(inplace=True)
    positions.to_csv(options.filtered_positions, sep="\t", header=False)


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
