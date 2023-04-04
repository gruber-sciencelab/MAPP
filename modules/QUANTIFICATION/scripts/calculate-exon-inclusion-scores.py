"""
##############################################################################
#
#   Calculate exon inclusion scores (forall exons, forall samples).
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 05-04-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
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
        "--transcript-quantification",
        dest="quant",
        required=True,
        help="Path to the merged transcript-level quantification results.",
    )
    parser.add_argument(
        "--exon-table",
        dest="exon_table",
        required=True,
        help="Path to the exon table with transcript-compatibility information.",
    )
    parser.add_argument(
        "--inclusion-table",
        dest="inclusion_table",
        required=True,
        help="Path for the output file with exon inclusion scores.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    # read the input tables
    expression = pd.read_csv(options.quant, sep="\t")
    columns = ["gene", "chr", "start", "end", "strand", "included", "total", "cluster"]
    exons = pd.read_csv(options.exon_table, sep="\t", names=columns)

    # create the final df: for every event calculate expression: included and total
    samples = list(expression.columns.values)
    columns_included = [s + "_included" for s in samples]
    columns_total = [s + "_total" for s in samples]
    columns = ["coordinates"] + columns_included + columns_total
    final_df = pd.DataFrame(columns=columns)
    for exon_id, row in exons.iterrows():
        included_sum = []  # sum of the expression of transcripts that include this exon
        for sample in samples:
            expr = sum(expression[sample].loc[row.included.split(",")])
            included_sum.append(expr)
        compatible_sum = (
            []
        )  # sum of the expression of transcripts that are compatible with this exon
        for sample in samples:
            expr = sum(expression[sample].loc[row.total.split(",")])
            compatible_sum.append(expr)
        final_df.loc[len(final_df)] = [exon_id] + included_sum + compatible_sum

    # save the inclusion table
    final_df.to_csv(options.inclusion_table, sep="\t", index=False)


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
