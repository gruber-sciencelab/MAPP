"""
##############################################################################
#
#   Calculate motif's activity differences Zscores
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 20-01-2020
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
import numpy as np


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
        "--activity-table",
        dest="activity_table",
        required=True,
        help="Path to the table with motifs activities and their stds.",
    )
    parser.add_argument(
        "--design-table",
        dest="design_table",
        required=True,
        help="Path to the design table.",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        required=True,
        help="Path for the output table with Z-scores.",
    )
    return parser


##############################################################################


def calculate_Zscores(A_b, design_table):
    """
    Calculate Zscores as the ratio: Act.Diff / Act.Diff.Std.
    """

    cols_list = []
    for s in design_table.index.values:
        cols_list.append("A_" + s)
    for s in design_table.index.values:
        cols_list.append("stdA_" + s)

    Zscores_df = pd.DataFrame(index=A_b.index.values)
    for c in cols_list:
        Zscores_df[c] = A_b[c]

    # calculate per-sample activity Z-scores
    for s in design_table.index.values:
        Zscores_df["zscores_" + s] = Zscores_df["A_" + s] / Zscores_df["stdA_" + s]
        cols_list.append("zscores_" + s)

    # calcualte motif Z-scores
    N_samples = len(design_table.index.values)
    for s in design_table.index.values:
        Zscores_df["Z^2_" + s] = Zscores_df["zscores_" + s] * Zscores_df["zscores_" + s]
    Zscores_df["sum_Z^2"] = 0.0
    for s in design_table.index.values:
        Zscores_df["sum_Z^2"] += Zscores_df["Z^2_" + s]
    Zscores_df["avg_Z^2"] = Zscores_df["sum_Z^2"] / N_samples
    Zscores_df["combined.Zscore"] = np.sqrt(Zscores_df["avg_Z^2"])

    # rename the columns
    new_cols_list = []
    for c in cols_list:
        c_prefix = c.split("_")[0]
        c_suffix = "_".join(c.split("_")[1:])
        if c_prefix == "A":
            new_cols_list.append("activities_" + c_suffix)
        if c_prefix == "stdA":
            new_cols_list.append("deltas_" + c_suffix)
        if c_prefix == "zscores":
            new_cols_list.append("zscores_" + c_suffix)

    # select columns for final table:
    cols_list.append("combined.Zscore")
    new_cols_list.append("combined.Zscore")
    Zscores_df = Zscores_df[cols_list].copy()
    Zscores_df.columns = new_cols_list
    return Zscores_df


def main():
    """Main body of the script."""

    # read the design table
    design_table = pd.read_csv(options.design_table, sep="\t", index_col=0)

    # read the table with activities and their std
    A_b_table = pd.read_csv(options.activity_table, sep="\t")
    A_b_table_fraction = A_b_table["fraction"].copy()
    del A_b_table["fraction"]

    # select motifs for which the model was not run
    not_fitted = A_b_table.loc[(A_b_table == 0).all(axis=1)].index.values

    # select a subtable containing motifs with activities fitted
    A_b_table = A_b_table.loc[~(A_b_table == 0).all(axis=1)]

    # calculate Z-scores
    Zscores_df = calculate_Zscores(A_b_table, design_table)

    # add the rows for motifs that was not fitted
    not_fitted_df = pd.DataFrame(
        np.nan, index=not_fitted, columns=Zscores_df.columns.values
    )
    Zscores_df = pd.concat([Zscores_df, not_fitted_df], axis=0)

    # save the table
    Zscores_df[
        "regions_present_fraction"
    ] = A_b_table_fraction  # restore the motif fractions
    Zscores_df.index.name = "ID"
    Zscores_df.to_csv(options.outfile, sep="\t", na_rep="NA")


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
