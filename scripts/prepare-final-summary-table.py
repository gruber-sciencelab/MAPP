"""
##############################################################################
#
#   Generate table with motifs' scores for the final HTML report.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 03-06-2021
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
import yaml
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
        "--ranking-score",
        dest="ranking_score",
        required=True,
        help="Motif ranking strategy (Max, Avg or RankProduct).",
    )
    parser.add_argument(
        "--summary-directory",
        dest="summary_directory",
        required=True,
        help="Path to the directory with MAPP summary.",
    )
    return parser


##############################################################################


def select_avgabs_score(ss3, ss5, pas):
    """Calculate Max-Abs Zscore over all sites"""
    if np.isnan(ss3):
        ss3 = 0
    else:
        ss3 = abs(ss3)
    if np.isnan(ss5):
        ss5 = 0
    else:
        ss5 = abs(ss5)

    if np.isnan(pas):
        pas = 0
    else:
        pas = abs(pas)
    return round((ss3 + ss5 + pas) / 3, 3)


def select_maxabs_score(scores):
    """Select the Max-Abs Zscore over all sites"""
    scores = [abs(i) for i in scores if not np.isnan(i)]
    if scores:
        return max(scores)
    else:
        return np.nan


def create_RankProduct_column(df):
    """Calculate per-motif Rank-Product over all three sites."""
    # Remark: all NaN values are ranked higher then actual scores but with the same rank(!)
    df["3'ss Zscore rank"] = df["3'ss Zscore"].rank(ascending=False, na_option="bottom")
    df["5'ss Zscore rank"] = df["5'ss Zscore"].rank(ascending=False, na_option="bottom")
    df["pas Zscore rank"] = df["pas Zscore"].rank(ascending=False, na_option="bottom")
    df["RankProduct"] = np.cbrt(
        df["3'ss Zscore rank"] * df["5'ss Zscore rank"] * df["pas Zscore rank"]
    )
    return df["RankProduct"]


def main():
    """Main body of the script."""

    # read in MAPP config file
    with open(os.path.join(options.summary_directory, "pipeline-config.yml")) as f:
        configfile = yaml.safe_load(f)
    if configfile["RES_sorting_strategy"] == "avg":
        ranking_column = "combined.standard.Zscore"
    else:
        assert configfile["RES_sorting_strategy"] == "max"
        ranking_column = "max.abs.standard.Zscore"

    # read in MAPP full result tables:
    ss3_table_original = pd.read_csv(
        os.path.join(options.summary_directory, "results", "results-3ss.tsv"),
        sep="\t",
        index_col=0,
    )
    ss5_table_original = pd.read_csv(
        os.path.join(options.summary_directory, "results", "results-5ss.tsv"),
        sep="\t",
        index_col=0,
    )
    pas_table_original = pd.read_csv(
        os.path.join(options.summary_directory, "results", "results-pas.tsv"),
        sep="\t",
        index_col=0,
    )

    # read in lists with selected motifs
    with open(
        os.path.join(options.summary_directory, "results", "motifs-list-3ss.txt")
    ) as f:
        ss3_motifs_list = f.read().splitlines()
    with open(
        os.path.join(options.summary_directory, "results", "motifs-list-5ss.txt")
    ) as f:
        ss5_motifs_list = f.read().splitlines()
    with open(
        os.path.join(options.summary_directory, "results", "motifs-list-pas.txt")
    ) as f:
        pas_motifs_list = f.read().splitlines()

    all_sign_motifs = list(
        set(ss3_motifs_list) | set(ss5_motifs_list) | set(pas_motifs_list)
    )

    # Subset only rows which denote motifs with SOME stat. sign. regulation:
    # Important assumption here: all three tables must have the same set of motifs
    # (that is why it was imprtant not to pre-filter)
    ss3_table = ss3_table_original[
        ss3_table_original["motif"].isin(all_sign_motifs)
    ].copy()
    ss5_table = ss5_table_original[
        ss5_table_original["motif"].isin(all_sign_motifs)
    ].copy()
    pas_table = pas_table_original[
        pas_table_original["motif"].isin(all_sign_motifs)
    ].copy()

    # subset only the columns we will need:
    cols_old = ["motif", ranking_column]
    ss3_table = ss3_table[cols_old].copy()
    ss3_table.columns = ["motif", "3'ss Zscore"]
    ss5_table = ss5_table[cols_old].copy()
    ss5_table.columns = ["motif", "5'ss Zscore"]
    pas_table = pas_table[cols_old].copy()
    pas_table.columns = ["motif", "pas Zscore"]

    # These tables are still sorted according to a respective sorting key.
    # I will select the "most significant" (max score) window by dropping duplicate motif-rows,
    # keeping only the FIRST occurance (iterating from the top)
    # ASSUMPTION: the highest-scoring window (max) is the one from which we capture statistical significance
    # This assumption holds since all these scores are already bg-normalized for each window,
    # therefore we can directly compare them.
    ss3_table.drop_duplicates(subset=["motif"], keep="first", inplace=True)
    ss5_table.drop_duplicates(subset=["motif"], keep="first", inplace=True)
    pas_table.drop_duplicates(subset=["motif"], keep="first", inplace=True)
    ss3_table.set_index("motif", drop=True, inplace=True)
    ss5_table.set_index("motif", drop=True, inplace=True)
    pas_table.set_index("motif", drop=True, inplace=True)

    # outer-join the three dataframes,
    # the new df is sorted alphabetically by the index name (motif)
    combined_df = pd.concat([ss3_table, ss5_table, pas_table], axis=1)
    combined_df.index.name = "motif"

    # at least one stat. sign. motif found
    if len(combined_df) > 0:

        # create a new ranking column:
        if options.ranking_score == "Max":
            combined_df["Ranking score"] = combined_df.apply(
                lambda row: select_maxabs_score(
                    [row["3'ss Zscore"], row["5'ss Zscore"], row["pas Zscore"]]
                ),
                axis=1,
            )
            combined_df.sort_values(by="Ranking score", inplace=True, ascending=False)
        elif options.ranking_score == "Avg":
            combined_df["Ranking score"] = combined_df.apply(
                lambda row: select_avgabs_score(
                    row["3'ss Zscore"], row["5'ss Zscore"], row["pas Zscore"]
                ),
                axis=1,
            )
            combined_df.sort_values(by="Ranking score", inplace=True, ascending=False)
        else:
            assert options.ranking_score == "RankProduct"
            combined_df["Ranking score"] = create_RankProduct_column(combined_df)
            combined_df.sort_values(by="Ranking score", inplace=True, ascending=True)

        # reorder columns:
        new_col_order = ["Ranking score", "3'ss Zscore", "5'ss Zscore", "pas Zscore"]
        combined_df = combined_df[new_col_order].copy()

        # add * marker near of motifs considered as stat. sign. for a given process
        combined_df = combined_df.round(decimals=3).applymap(str)
        for m, row in combined_df.iterrows():
            if m in ss3_motifs_list:
                combined_df.at[m, "3'ss Zscore"] = (
                    combined_df.at[m, "3'ss Zscore"] + " *"
                )
            if m in ss5_motifs_list:
                combined_df.at[m, "5'ss Zscore"] = (
                    combined_df.at[m, "5'ss Zscore"] + " *"
                )
            if m in pas_motifs_list:
                combined_df.at[m, "pas Zscore"] = combined_df.at[m, "pas Zscore"] + " *"

        # add a column with relative paths to the PNG heatmaps:
        combined_df["Activity Map"] = combined_df.apply(
            lambda row: os.path.join(
                ".", "results", "heatmaps", row.name, "heatmap.png"
            ),
            axis=1,
        )

        # if this was a run in pwm mode - add relative paths to the seqlogos:
        if configfile["CSM_matrix_type"] == "pwms":
            combined_df["motif"] = combined_df.index.values
            new_col_order = [
                "motif",
                "Ranking score",
                "3'ss Zscore",
                "5'ss Zscore",
                "pas Zscore",
                "Activity Map",
            ]
            combined_df = combined_df[new_col_order].copy()
            combined_df["Sequence Logos"] = combined_df.apply(
                lambda row: os.path.join(
                    ".", "sequence-logos", row.name + ".png"
                ),
                axis=1,
            )
            combined_df = combined_df.set_index("Sequence Logos")

        # save the final table:
        combined_df.to_csv(
            os.path.join(options.summary_directory, "main-table.tsv"),
            sep="\t",
            na_rep="NA",
        )

    else:  # edge case if no stat. sign. motifs found at all
        new_col_order = [
            "motif",
            "Ranking score",
            "3'ss Zscore",
            "5'ss Zscore",
            "pas Zscore",
            "Activity Map",
        ]
        if configfile["CSM_matrix_type"] == "pwms":
            new_col_order = ["Sequence Logos"] + new_col_order
        combined_df = pd.DataFrame(columns=new_col_order)
        # save the final table:
        combined_df.to_csv(
            os.path.join(options.summary_directory, "main-table.tsv"),
            sep="\t",
            na_rep="NA",
            index=False,
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
