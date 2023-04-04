"""
##############################################################################
#
#   Plot heatmaps of samples' Z-scores in every window over every site
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
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns

mpl.use("Agg")


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
        "--tables-main-directory",
        dest="indir",
        required=True,
        help="Path to the input directory with Z-scores tables per motif.",
    )
    parser.add_argument(
        "--output-directory",
        dest="outdir",
        required=True,
        help="Path for the output directory with Z-scores plots per motif.",
    )
    parser.add_argument(
        "--dpi",
        dest="dpi",
        default=300,
        help="Dots per inch parameter for the PNG plots (default: 300).",
    )
    return parser


##############################################################################

# GLOBAL SETTINGS FOR THE PLOTS:
PLOT_figsize_height = 20
PLOT_figsize_width = 8
PLOT_width_ratios = (10, 10, 10, 1)
PLOT_height_ratios = (7, 1)
PLOT_right_colorbar_label_fontsize = 20
PLOT_right_colorbar_tick_width = 2
PLOT_subtitles_fontsize = 20
PLOT_subtitles_padding = 40

##############################################################################


def plot_heatmap(m, m_path):
    """Build the heatmap Z-scores summary plot for a given motif."""

    # read the table with standardized Z-scores
    path_3ss = os.path.join(m_path, "Zscores_3ss.tsv")
    path_5ss = os.path.join(m_path, "Zscores_5ss.tsv")
    path_pas = os.path.join(m_path, "Zscores_pas.tsv")
    m_df_3ss = pd.read_csv(path_3ss, sep="\t", index_col=0)
    m_df_5ss = pd.read_csv(path_5ss, sep="\t", index_col=0)
    m_df_pas = pd.read_csv(path_pas, sep="\t", index_col=0)

    # read the file with statistically significant columns
    path_3ss = os.path.join(m_path, "statistically_significant_columns_3ss.tsv")
    path_5ss = os.path.join(m_path, "statistically_significant_columns_5ss.tsv")
    path_pas = os.path.join(m_path, "statistically_significant_columns_pas.tsv")
    with open(path_3ss) as f:
        m_columns_3ss = f.read().splitlines()
    with open(path_5ss) as f:
        m_columns_5ss = f.read().splitlines()
    with open(path_pas) as f:
        m_columns_pas = f.read().splitlines()

    # read the file with statistically significant cells
    path_3ss = os.path.join(m_path, "statistically_significant_cells_3ss.tsv")
    path_5ss = os.path.join(m_path, "statistically_significant_cells_5ss.tsv")
    path_pas = os.path.join(m_path, "statistically_significant_cells_pas.tsv")
    with open(path_3ss) as f:
        m_cells_3ss = [i.replace("\t", "$") for i in f.read().splitlines()]
    with open(path_5ss) as f:
        m_cells_5ss = [i.replace("\t", "$") for i in f.read().splitlines()]
    with open(path_pas) as f:
        m_cells_pas = [i.replace("\t", "$") for i in f.read().splitlines()]

    # construct the df's with stat. sign. annotation for the heatmaps
    annot_df_3ss = m_df_3ss.copy().astype(str)
    for row in annot_df_3ss.index.values:
        for col in annot_df_3ss.columns.values:
            if col + "$" + row in m_cells_3ss:
                annot_df_3ss.at[row, col] = "*"
            else:
                annot_df_3ss.at[row, col] = ""
    annot_df_5ss = m_df_5ss.copy().astype(str)
    for row in annot_df_5ss.index.values:
        for col in annot_df_5ss.columns.values:
            if col + "$" + row in m_cells_5ss:
                annot_df_5ss.at[row, col] = "*"
            else:
                annot_df_5ss.at[row, col] = ""
    annot_df_pas = m_df_pas.copy().astype(str)
    for row in annot_df_pas.index.values:
        for col in annot_df_pas.columns.values:
            if col + "$" + row in m_cells_pas:
                annot_df_pas.at[row, col] = "*"
            else:
                annot_df_pas.at[row, col] = ""

    # adjust the size of text elements based on samples/windows
    _3ss_n_samples = len(m_df_3ss.index.values)
    _5ss_n_samples = len(m_df_5ss.index.values)
    _pas_n_samples = len(m_df_pas.index.values)
    assert _3ss_n_samples == _5ss_n_samples == _pas_n_samples
    _n_samples = _pas_n_samples
    _3ss_n_windows = len(m_df_3ss.columns.values)
    _5ss_n_windows = len(m_df_5ss.columns.values)
    _pas_n_windows = len(m_df_pas.columns.values)
    max_n_windows = max(_3ss_n_windows, _5ss_n_windows, _pas_n_windows)
    max_samples_windows = max(_n_samples, max_n_windows)
    x_ticks_size = 5.0 if max_n_windows >= 50 else 10 - np.floor(max_n_windows / 10)
    y_ticks_size = 0.5 if _n_samples >= 100 else 10 - np.floor(_n_samples / 10)
    annot_star_size = (
        0.5 if max_samples_windows >= 100 else 10 - np.floor(max_samples_windows / 10)
    )

    # create the axes layout
    fig, ax = mpl.pyplot.subplots(
        nrows=2,
        ncols=4,
        figsize=(PLOT_figsize_height, PLOT_figsize_width),
        gridspec_kw={
            "width_ratios": PLOT_width_ratios,
            "height_ratios": PLOT_height_ratios,
        },
    )

    # extract the max abs standard Z-score
    m_df_3ss_max = m_df_3ss.abs().max(skipna=True).max(skipna=True)
    m_df_5ss_max = m_df_5ss.abs().max(skipna=True).max(skipna=True)
    m_df_pas_max = m_df_pas.abs().max(skipna=True).max(skipna=True)
    max_value = pd.Series([m_df_3ss_max, m_df_5ss_max, m_df_pas_max]).dropna().max()

    # 3ss heatmap
    g_3ss = sns.heatmap(
        m_df_3ss,
        ax=ax[0, 0],
        cbar=False,
        annot=annot_df_3ss,
        fmt="",
        annot_kws={"size": annot_star_size, "va": "center_baseline"},
        cmap="coolwarm",
        linewidth=0.5,
        linecolor="black",
        vmin=-max_value,
        vmax=max_value,
        xticklabels=True,
        yticklabels=True,
    )
    g_3ss.set_facecolor("white")
    mpl.rcParams["hatch.linewidth"] = 0.5
    mpl.rcParams["hatch.color"] = "grey"
    ax[0, 0].pcolor(np.where(m_df_3ss.isna(), 0, np.nan), hatch="/", alpha=0)
    ax[0, 0].set_title(
        "3'ss", fontsize=PLOT_subtitles_fontsize, pad=PLOT_subtitles_padding
    )
    ax[0, 0].tick_params(
        axis="x", labelsize=x_ticks_size,
    )
    ax[0, 0].tick_params(
        axis="y", labelsize=y_ticks_size, labelrotation=1,
    )
    ax[0, 0].set_xticklabels(ax[0, 0].get_xticklabels(), rotation=90)
    # annotate statistical significancy for the columns:
    marker = ["*" if win in m_columns_3ss else "" for win in m_df_3ss.columns.values]
    ax_0_1_copy = ax[0, 0].twiny()
    ax_0_1_copy.set_xlim([0, ax[0, 0].get_xlim()[1]])
    ax_0_1_copy.set_xticks(ax[0, 0].get_xticks())
    ax_0_1_copy.set_xticklabels(marker, fontsize=2 + x_ticks_size)
    ax_0_1_copy.tick_params(top=False)

    # 5ss heatmap
    g_5ss = sns.heatmap(
        m_df_5ss,
        ax=ax[0, 1],
        cbar=False,
        annot=annot_df_5ss,
        fmt="",
        annot_kws={"size": annot_star_size, "va": "center_baseline"},
        cmap="coolwarm",
        linewidth=0.5,
        linecolor="black",
        vmin=-max_value,
        vmax=max_value,
        xticklabels=True,
        yticklabels=True,
    )
    g_5ss.set_facecolor("white")
    mpl.rcParams["hatch.linewidth"] = 0.5
    mpl.rcParams["hatch.color"] = "grey"
    ax[0, 1].pcolor(np.where(m_df_5ss.isna(), 0, np.nan), hatch="/", alpha=0)
    ax[0, 1].set_yticks([])
    ax[0, 1].set_title(
        "5'ss", fontsize=PLOT_subtitles_fontsize, pad=PLOT_subtitles_padding
    )
    ax[0, 1].tick_params(
        axis="x", labelsize=x_ticks_size,
    )
    ax[0, 1].set_xticklabels(ax[0, 1].get_xticklabels(), rotation=90)
    # annotate statistical significancy for the columns:
    marker = ["*" if win in m_columns_5ss else "" for win in m_df_5ss.columns.values]
    ax_0_2_copy = ax[0, 1].twiny()
    ax_0_2_copy.set_xlim([0, ax[0, 1].get_xlim()[1]])
    ax_0_2_copy.set_xticks(ax[0, 1].get_xticks())
    ax_0_2_copy.set_xticklabels(marker, fontsize=2 + x_ticks_size)
    ax_0_2_copy.tick_params(top=False)

    # pas heatmap
    g_pas = sns.heatmap(
        m_df_pas,
        ax=ax[0, 2],
        cbar=False,
        annot=annot_df_pas,
        fmt="",
        annot_kws={"size": annot_star_size, "va": "center_baseline"},
        cmap="coolwarm",
        linewidth=0.5,
        linecolor="black",
        vmin=-max_value,
        vmax=max_value,
        xticklabels=True,
        yticklabels=True,
    )
    g_pas.set_facecolor("white")
    mpl.rcParams["hatch.linewidth"] = 0.5
    mpl.rcParams["hatch.color"] = "grey"
    ax[0, 2].pcolor(np.where(m_df_pas.isna(), 0, np.nan), hatch="/", alpha=0)
    ax[0, 2].set_yticks([])
    ax[0, 2].set_title(
        "pas", fontsize=PLOT_subtitles_fontsize, pad=PLOT_subtitles_padding
    )
    ax[0, 2].tick_params(
        axis="x", labelsize=x_ticks_size,
    )
    ax[0, 2].set_xticklabels(ax[0, 2].get_xticklabels(), rotation=90)
    # annotate statistical significancy for the columns:
    marker = ["*" if win in m_columns_pas else "" for win in m_df_pas.columns.values]
    ax_0_3_copy = ax[0, 2].twiny()
    ax_0_3_copy.set_xlim([0, ax[0, 2].get_xlim()[1]])
    ax_0_3_copy.set_xticks(ax[0, 2].get_xticks())
    ax_0_3_copy.set_xticklabels(marker, fontsize=2 + x_ticks_size)
    ax_0_3_copy.tick_params(top=False)

    # right color bar
    mpl.pyplot.colorbar(
        mpl.pyplot.cm.ScalarMappable(
            cmap="coolwarm", norm=mpl.colors.Normalize(vmin=-max_value, vmax=max_value),
        ),
        cax=ax[0, 3],
    )
    ax[0, 3].yaxis.set_ticks_position("right")
    ax[0, 3].tick_params(
        axis="y",
        labelsize=PLOT_right_colorbar_label_fontsize,
        width=PLOT_right_colorbar_tick_width,
        labelrotation=1,
    )

    # clear the bottom axes
    ax[1, 0].axis("off")
    ax[1, 1].axis("off")
    ax[1, 2].axis("off")
    ax[1, 3].axis("off")

    # save the plot
    m_dir = os.path.join(options.outdir, m)
    os.mkdir(m_dir)
    mpl.pyplot.savefig(
        os.path.join(m_dir, "heatmap.png"), format="png", dpi=int(options.dpi)
    )
    mpl.pyplot.savefig(os.path.join(m_dir, "heatmap.pdf"), format="pdf")
    mpl.pyplot.close()


def main():
    """Main body of the script."""

    # create main output directory
    os.mkdir(options.outdir)

    # iterate & plot over every motif:
    for m_path in glob.glob(os.path.join(options.indir, "*/")):
        m = m_path.split("/")[-2]
        plot_heatmap(m, m_path)


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
