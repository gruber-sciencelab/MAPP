"""
##############################################################################
#
#   Prepare tables for samples' Z-scores heatmaps (all windows/all sites)
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
        "--samples-design-table",
        dest="design_table",
        required=True,
        help="Path to the TSV file with information about samples.",
    )
    parser.add_argument(
        "--max-pvalue",
        dest="maxpval",
        required=True,
        help="Maximal corrected p-value for statistical significancy.",
    )
    parser.add_argument(
        "--sorting-strategy",
        dest="sorting_strategy",
        required=True,
        help="Strategy utilised for sorting the results: avg OR max.",
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
        "--infile-splicing-3ss-motifs",
        dest="motifs_3ss",
        required=True,
        help="Path to the text file with top motifs (3ss).",
    )
    parser.add_argument(
        "--infile-splicing-5ss-motifs",
        dest="motifs_5ss",
        required=True,
        help="Path to the text file with top motifs (5ss).",
    )
    parser.add_argument(
        "--infile-polyadenylation-pas-motifs",
        dest="motifs_pas",
        required=True,
        help="Path to the text file with top motifs (pas).",
    )
    parser.add_argument(
        "--output-directory",
        dest="outdir",
        required=True,
        help="Path for the output directory with Z-scores tables.",
    )
    return parser


##############################################################################


def sort_windows(l):
    """Sort a list of sliding windows."""
    parsed_l = [w.split(".") for w in l]
    parsed_l = [
        (
            int(w[0][1:].split("to")[0]),
            int(w[0][1:].split("to")[1]),
            int(w[1][1:].split("to")[0]),
            int(w[1][1:].split("to")[1]),
        )
        for w in parsed_l
    ]
    parsed_l = sorted(parsed_l, key=lambda tup: (-tup[0], tup[3]))
    new_l = []
    for w in parsed_l:
        new_l.append(
            "u" + str(w[0]) + "to" + str(w[1]) + ".d" + str(w[2]) + "to" + str(w[3])
        )
    return new_l


def prepare_Zscores_table_per_motif_per_site(df, m, samples, maxpval):
    """Create a table with per-window-per-sample Z-scores for a given site."""
    m_df = df[df["motif"] == m].copy()
    m_df = m_df.set_index("region")

    # sort the windows and the whole df accordingly
    sorted_windows = sort_windows(m_df.index.values)
    m_df = m_df.loc[sorted_windows]

    # prepare df for the per-sample, per-window standardized Z-scores
    Zscores_table = pd.DataFrame(columns=sorted_windows)

    # select per-sample standardized Z-scores
    for s in samples:
        Zscores_table.loc[s] = m_df["standard_zscores_" + s]

    # =========================================================
    # select which windows/samples are statistically significant

    pval_cols = ["corr_pval_" + s for s in samples]

    sign_cells = []
    for row in m_df.index.values:
        for col in pval_cols:
            sample_ID = col.split("corr_pval_")[-1]
            if m_df.at[row, col] <= maxpval:
                sign_cells.append(row + "\t" + sample_ID)

    sign_columns = []
    for row in m_df.index.values:
        pvals = []
        for col in pval_cols:
            pvals.append(m_df.at[row, col])
        pvals = pd.Series(pvals)
        pvals = pvals <= maxpval
        if options.sorting_strategy == "avg":
            if pvals.all():
                sign_columns.append(row)
        else:
            assert options.sorting_strategy == "max"
            if pvals.any():
                sign_columns.append(row)

    return (Zscores_table, sign_cells, sign_columns)


def main():
    """Main body of the script."""

    # read the input data
    design_table = pd.read_csv(options.design_table, sep="\t", index_col=0)
    df_3ss = pd.read_csv(options.results_3ss, sep="\t", index_col=0)
    df_5ss = pd.read_csv(options.results_5ss, sep="\t", index_col=0)
    df_pas = pd.read_csv(options.results_pas, sep="\t", index_col=0)
    with open(options.motifs_3ss) as f:
        motifs_3ss = f.read().splitlines()
    with open(options.motifs_5ss) as f:
        motifs_5ss = f.read().splitlines()
    with open(options.motifs_pas) as f:
        motifs_pas = f.read().splitlines()

    # create the output directory
    os.mkdir(options.outdir)

    all_motifs = set(motifs_3ss + motifs_5ss + motifs_pas)
    for m in all_motifs:

        # create a directory for a given motif
        os.mkdir(os.path.join(options.outdir, m))

        # extract info from the tables
        (
            m_df_3ss,
            sign_cells_3ss,
            sign_columns_3ss,
        ) = prepare_Zscores_table_per_motif_per_site(
            df_3ss, m, design_table.index.values, float(options.maxpval)
        )
        (
            m_df_5ss,
            sign_cells_5ss,
            sign_columns_5ss,
        ) = prepare_Zscores_table_per_motif_per_site(
            df_5ss, m, design_table.index.values, float(options.maxpval)
        )
        (
            m_df_pas,
            sign_cells_pas,
            sign_columns_pas,
        ) = prepare_Zscores_table_per_motif_per_site(
            df_pas, m, design_table.index.values, float(options.maxpval)
        )

        # save standardized Z-scores DataFrames
        m_df_3ss.to_csv(
            os.path.join(options.outdir, m, "Zscores_3ss.tsv"),
            sep="\t",
            index_label=False,
        )
        m_df_5ss.to_csv(
            os.path.join(options.outdir, m, "Zscores_5ss.tsv"),
            sep="\t",
            index_label=False,
        )
        m_df_pas.to_csv(
            os.path.join(options.outdir, m, "Zscores_pas.tsv"),
            sep="\t",
            index_label=False,
        )

        # save TSV files with statistically significant cells for the heatmap
        with open(
            os.path.join(options.outdir, m, "statistically_significant_cells_3ss.tsv"),
            "w",
        ) as f:
            for i in sign_cells_3ss:
                f.write(i + os.linesep)
        with open(
            os.path.join(options.outdir, m, "statistically_significant_cells_5ss.tsv"),
            "w",
        ) as f:
            for i in sign_cells_5ss:
                f.write(i + os.linesep)
        with open(
            os.path.join(options.outdir, m, "statistically_significant_cells_pas.tsv"),
            "w",
        ) as f:
            for i in sign_cells_pas:
                f.write(i + os.linesep)

        # save TXT files with statistically significant columns for the heatmap
        with open(
            os.path.join(
                options.outdir, m, "statistically_significant_columns_3ss.tsv"
            ),
            "w",
        ) as f:
            for i in sign_columns_3ss:
                f.write(i + os.linesep)
        with open(
            os.path.join(
                options.outdir, m, "statistically_significant_columns_5ss.tsv"
            ),
            "w",
        ) as f:
            for i in sign_columns_5ss:
                f.write(i + os.linesep)
        with open(
            os.path.join(
                options.outdir, m, "statistically_significant_columns_pas.tsv"
            ),
            "w",
        ) as f:
            for i in sign_columns_pas:
                f.write(i + os.linesep)


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
