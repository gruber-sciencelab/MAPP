"""
##############################################################################
#
#   Normalize pas expression for per-sample library size (result in TPMs)
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
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
        "--tandem-pas-expression",
        dest="tandem_pas_expression",
        required=True,
        help="Path to the TSV file with quantified tandem pas expression values.",
    )
    parser.add_argument(
        "--distal-pas-expression",
        dest="distal_pas_expression",
        required=True,
        help="Path to the TSV file with quantified distal pas expression values.",
    )
    parser.add_argument(
        "--expression-pseudocount",
        dest="expression-pseudocount",
        required=True,
        help="Pseudocount constant to add to the expression table.",
    )
    parser.add_argument(
        "--normalized-tandem-pas-expression-original",
        dest="normalized_tandem_pas_expression_original",
        required=True,
        help="Path for the TSV file with normalized tandem pas expression values (original).",
    )
    parser.add_argument(
        "--normalized-tandem-pas-expression-pseudocount",
        dest="normalized_tandem_pas_expression_pseudocount",
        required=True,
        help="Path for the TSV file with normalized tandem pas expression values (pseudocount).",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    tandem_pas_expression = pd.read_csv(
        options.tandem_pas_expression, sep="\t"
    ).replace(to_replace=-1.0, value=np.nan)

    distal_pas_expression = pd.read_csv(
        options.distal_pas_expression, sep="\t"
    ).replace(to_replace=-1.0, value=np.nan)

    samples = tandem_pas_expression.columns.values[10:]
    total_expression = {s: None for s in samples}

    for s in samples:
        # Watch out for internal type conversion, we don't want to loose precision!
        # https://stackoverflow.com/a/46979587/2340598
        total_expression[s] = tandem_pas_expression[s].astype(np.float64).sum(
            skipna=True
        ) + distal_pas_expression[s].astype(np.float64).sum(skipna=True)
        tandem_pas_expression[s] = (
            tandem_pas_expression[s].astype(np.float64) / total_expression[s] * 1000000
        )

    tandem_pas_expression.replace(to_replace=np.nan, value=-1.0).to_csv(
        options.normalized_tandem_pas_expression_original,
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    tandem_pas_expression.replace(to_replace=np.nan, value=0.0, inplace=True)
    for s in samples:
        tandem_pas_expression[s] = tandem_pas_expression[s] + 1.0
    tandem_pas_expression.to_csv(
        options.normalized_tandem_pas_expression_pseudocount,
        sep="\t",
        index=False,
        float_format="%.6f",
    )


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
