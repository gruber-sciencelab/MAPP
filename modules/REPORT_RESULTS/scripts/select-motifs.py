"""
##############################################################################
#
#   Select top N distinct motifs with highest (statistically significant)
#   activity Z-score (for every site separately)
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 04-06-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import os
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
        "--topN-motifs",
        dest="N",
        default=1000000,  # by default: effectively select all stat. sign. motifs
        required=False,
        help="Number of top motifs to select.",
    )
    parser.add_argument(
        "--infile-splicing-3ss",
        dest="results_3ss",
        required=True,
        help="Annotated results table (3ss).",
    )
    parser.add_argument(
        "--infile-splicing-5ss",
        dest="results_5ss",
        required=True,
        help="Annotated results table (5ss).",
    )
    parser.add_argument(
        "--infile-polyadenylation-pas",
        dest="results_pas",
        required=True,
        help="Annotated results table (pas).",
    )
    parser.add_argument(
        "--outfile-splicing-3ss-motifs",
        dest="motifs_3ss",
        required=True,
        help="Path for the text file with top motifs (3ss).",
    )
    parser.add_argument(
        "--outfile-splicing-5ss-motifs",
        dest="motifs_5ss",
        required=True,
        help="Path for the text file with top motifs (5ss).",
    )
    parser.add_argument(
        "--outfile-polyadenylation-pas-motifs",
        dest="motifs_pas",
        required=True,
        help="Path for the text file with top motifs (pas).",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    df = pd.read_csv(options.results_3ss, sep="\t", index_col=0)
    df = df[df["significance-marker"]]
    motifs = []
    for ID, row in df.iterrows():
        if len(motifs) == int(options.N):
            break
        m = ID.split("|")[-1]
        if m not in motifs:
            motifs.append(m)
    with open(options.motifs_3ss, "w") as f:
        for m in motifs:
            f.write(m + os.linesep)

    df = pd.read_csv(options.results_5ss, sep="\t", index_col=0)
    df = df[df["significance-marker"]]
    motifs = []
    for ID, row in df.iterrows():
        if len(motifs) == int(options.N):
            break
        m = ID.split("|")[-1]
        if m not in motifs:
            motifs.append(m)
    with open(options.motifs_5ss, "w") as f:
        for m in motifs:
            f.write(m + os.linesep)

    df = pd.read_csv(options.results_pas, sep="\t", index_col=0)
    df = df[df["significance-marker"]]
    motifs = []
    for ID, row in df.iterrows():
        if len(motifs) == int(options.N):
            break
        m = ID.split("|")[-1]
        if m not in motifs:
            motifs.append(m)
    with open(options.motifs_pas, "w") as f:
        for m in motifs:
            f.write(m + os.linesep)


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
