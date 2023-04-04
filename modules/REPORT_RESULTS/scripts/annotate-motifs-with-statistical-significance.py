"""
##############################################################################
#
#   Annotate statistical significancy for all motifs in all windows
#   around all sites based on the provided cutoff parameter.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
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
        "--sorting-strategy",
        dest="sorting_strategy",
        required=True,
        help="Strategy utilised for sorting the results: avg OR max.",
    )
    parser.add_argument(
        "--max-pvalue",
        dest="maxpval",
        required=True,
        help="Maximal corrected p-value for statistical significancy.",
    )
    parser.add_argument(
        "--infile-splicing-3ss",
        dest="results_3ss",
        required=True,
        help="Results table (3ss).",
    )
    parser.add_argument(
        "--infile-splicing-5ss",
        dest="results_5ss",
        required=True,
        help="Results table (5ss).",
    )
    parser.add_argument(
        "--infile-polyadenylation-pas",
        dest="results_pas",
        required=True,
        help="Results table (pas).",
    )
    parser.add_argument(
        "--outfile-splicing-3ss",
        dest="results_3ss_updated",
        required=True,
        help="Path for the TSV table with annotated information (sorted).",
    )
    parser.add_argument(
        "--outfile-splicing-5ss",
        dest="results_5ss_updated",
        required=True,
        help="Path for the TSV table with annotated information (sorted).",
    )
    parser.add_argument(
        "--outfile-polyadenylation-pas",
        dest="results_pas_updated",
        required=True,
        help="Path for the TSV table with annotated information (sorted).",
    )
    return parser


##############################################################################


def avg_stat_sign(row):
    """
    In the "avg" mode a motif is stat. sign. in a given window iff
    the pval is below a cutoff in AT LEAST HALF of the samples.
    """
    pval_prefix = "corr_pval"
    pval_names = [i for i in list(row.keys()) if i.startswith(pval_prefix)]
    below_cutoff_counter = 0
    for pname in pval_names:
        if np.isnan(row[pname]):
            return False
        if row[pname] < float(options.maxpval):
            below_cutoff_counter += 1
    return below_cutoff_counter >= len(pval_names) / 2


def max_stat_sign(row):
    """
    In the "max" mode a motif is stat. sign. in a given window iff
    the pval from the max-abs-Zscore is below the cutoff.
    """
    pval_prefix = "corr_pval"
    pval_names = [i for i in list(row.keys()) if i.startswith(pval_prefix)]
    pvals = [row[pname] for pname in pval_names]
    for pval in pvals:
        if np.isnan(pval):
            return False
    # max abs Zscore corresponds to the min pval
    minpval = min(pvals)
    return minpval <= float(options.maxpval)


def main():
    """Main body of the script."""

    # In the "avg" mode a motif is stat. sign. in a given window iff
    # the pval is below a cutoff in ALL samples.
    if options.sorting_strategy == "avg":

        df = pd.read_csv(options.results_3ss, sep="\t", index_col=0)
        df["significance-marker"] = df.apply(lambda row: avg_stat_sign(row), axis=1)
        df.to_csv(options.results_3ss_updated, sep="\t")

        df = pd.read_csv(options.results_5ss, sep="\t", index_col=0)
        df["significance-marker"] = df.apply(lambda row: avg_stat_sign(row), axis=1)
        df.to_csv(options.results_5ss_updated, sep="\t")

        df = pd.read_csv(options.results_pas, sep="\t", index_col=0)
        df["significance-marker"] = df.apply(lambda row: avg_stat_sign(row), axis=1)
        df.index.name = "ID"  # inconsistent output between KAPACv2 and MAEI
        # correct inconsistent column order between KAPACv2 and MAEI
        colmarker = list(df.columns.values).index("regions_present_fraction")
        newcolorder = (
            list(df.columns.values[1:colmarker])
            + [df.columns.values[0]]
            + list(df.columns.values[colmarker:])
        )
        df = df[newcolorder]
        df.to_csv(options.results_pas_updated, sep="\t")

    # In the "max" mode a motif is stat. sign. in a given window iff
    # the pval from the max-abs-Zscore is below the cutoff.
    else:
        assert options.sorting_strategy == "max"

        df = pd.read_csv(options.results_3ss, sep="\t", index_col=0)
        df["significance-marker"] = df.apply(lambda row: max_stat_sign(row), axis=1)
        df.to_csv(options.results_3ss_updated, sep="\t")

        df = pd.read_csv(options.results_5ss, sep="\t", index_col=0)
        df["significance-marker"] = df.apply(lambda row: max_stat_sign(row), axis=1)
        df.to_csv(options.results_5ss_updated, sep="\t")

        df = pd.read_csv(options.results_pas, sep="\t", index_col=0)
        df["significance-marker"] = df.apply(lambda row: max_stat_sign(row), axis=1)
        df.index.name = "ID"  # inconsistent output between KAPACv2 and MAEI
        # correct inconsistent column order between KAPACv2 and MAEI
        colmarker = list(df.columns.values).index("regions_present_fraction")
        newcolorder = (
            list(df.columns.values[1:colmarker])
            + [df.columns.values[0]]
            + list(df.columns.values[colmarker:])
        )
        df = df[newcolorder]
        df.to_csv(options.results_pas_updated, sep="\t")


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
