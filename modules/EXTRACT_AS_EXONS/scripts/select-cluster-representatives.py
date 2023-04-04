"""
##############################################################################
#
#   Select exon representative subset according to n.o.
#   transcripts that contain them.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 26-03-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
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
        "--infile",
        dest="infile",
        required=True,
        help="Path to the file with clustered exons.",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        required=True,
        help="Path for the file with representative exons.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    clustered = pd.read_csv(options.infile, sep="\t", header=None)
    clustered.columns = [
        "exon",
        "gene",
        "chr",
        "start",
        "end",
        "strand",
        "T_including",
        "T_total",
        "cluster",
    ]
    clustered["no_including_transcripts"] = clustered.apply(
        lambda x: len(x["T_including"].split(",")), axis=1
    )

    # select exons supported by the highest number of transcripts
    dominant_exons = []
    for i, df in clustered.groupby(by="cluster"):
        exon = df.at[int(df[["no_including_transcripts"]].idxmax()), "exon"]
        dominant_exons.append(exon)
    filtered_df = clustered[clustered["exon"].isin(dominant_exons)].copy()
    del filtered_df["no_including_transcripts"]
    filtered_df.to_csv(options.outfile, sep="\t", header=False, index=False)


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
