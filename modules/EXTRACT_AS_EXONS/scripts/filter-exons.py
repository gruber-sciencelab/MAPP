"""
##############################################################################
#
#   Filter the protein-coding set of exons that have valid(long-enough)
#   regions around 5/3-splice-site.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 26-03-2020
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import os
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
from GTF import *
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
        "--in",
        dest="input",
        required=True,
        help="Path to the SUPPA file with SE events.",
    )
    parser.add_argument(
        "--gtf",
        dest="gtf",
        required=True,
        help="Path to the annotation file in gtf format.",
    )
    parser.add_argument(
        "--biotype_filters",
        dest="biotype_filters",
        required=True,
        help="Comma-separated list of biotypes to filter the gtf ('all' string turns off the filtering).",
    )
    parser.add_argument(
        "--out", dest="output", required=True, help="Path for the output file."
    )
    parser.add_argument(
        "--region_size_3ss_up",
        dest="region_size_3ss_up",
        required=True,
        help="Size of the region upstream 3' splicing site.",
    )
    parser.add_argument(
        "--region_size_3ss_down",
        dest="region_size_3ss_down",
        required=True,
        help="Size of the region downstream 3' splicing site.",
    )
    parser.add_argument(
        "--region_size_5ss_up",
        dest="region_size_5ss_up",
        required=True,
        help="Size of the region upstream 5' splicing site.",
    )
    parser.add_argument(
        "--region_size_5ss_down",
        dest="region_size_5ss_down",
        required=True,
        help="Size of the region downstream 5' splicing site.",
    )
    parser.add_argument(
        "--min_exon_len",
        dest="min_exon_len",
        required=True,
        help="Minimal exon length filter.",
    )
    parser.add_argument(
        "--max_exon_len",
        dest="max_exon_len",
        required=True,
        help="Maximal exon length filter.",
    )
    parser.add_argument(
        "--fai",
        dest="fai",
        required=True,
        help="Path to the genome index in fai format.",
    )
    parser.add_argument(
        "--filtering-logfile",
        dest="filter_log",
        required=True,
        help="Path to the file for filtering logs.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    # get dict of chrom. lengths
    with open(options.fai, "r") as fai:
        chrom_len_dict = {line.split("\t")[0]: int(line.split("\t")[1]) for line in fai}

    df_gtf = dataframe(options.gtf)  # from GTF.py
    df_exons = pd.read_csv(options.input, sep="\t")
    df_exons["ID"] = df_exons.apply(
        lambda x: str(x.seqname)
        + ":"
        + x.event_id.split("-")[1].replace(":", "-")
        + ":"
        + x.event_id[-1],
        axis=1,
    )

    with open(options.filter_log, "w") as log:

        log.write("filter-exons.py" + os.linesep)
        log.write(str(df_gtf.shape[0]) + "\tentries in GTF file" + os.linesep)

        # get rid of duplicates (that differ on other fields) here.
        if options.biotype_filters != "all":
            df_gtf = (
                df_gtf[["seqname", "start", "end", "strand", "feature", "gene_biotype"]]
                .copy()
                .drop_duplicates()
            )
        else:
            df_gtf = (
                df_gtf[["seqname", "start", "end", "strand", "feature"]]
                .copy()
                .drop_duplicates()
            )
        df_gtf["name"] = (
            df_gtf["seqname"]
            + ":"
            + df_gtf["start"]
            + "-"
            + df_gtf["end"]
            + ":"
            + df_gtf["strand"]
        )

        # get only exons and filter on gene_biotype (optionally)
        df_gtf = df_gtf[df_gtf["feature"] == "exon"].copy()
        log.write(str(df_gtf.shape[0]) + "\tunique exons in GTF" + os.linesep)
        if options.biotype_filters != "all":
            biotypes = options.biotype_filters.split(",")
            df_gtf = df_gtf[df_gtf["gene_biotype"].isin(biotypes)].copy()
            log.write(
                str(df_gtf.shape[0])
                + "\texons filtered on the gene biotype"
                + os.linesep
            )
        else:
            log.write("gene_biotype filtering turned off" + os.linesep)
        df_gtf[["start", "end"]] = df_gtf[["start", "end"]].astype(int)

        # filter on minimal length
        df_gtf["len"] = df_gtf["end"] - df_gtf["start"]
        df_gtf = df_gtf[
            df_gtf["len"]
            > int(options.min_exon_len)
        ].copy()
        log.write(
            str(df_gtf.shape[0])
            + "\tlength > min length parameter"
            + os.linesep
        )

        # filter on maximal length
        df_gtf["len"] = df_gtf["end"] - df_gtf["start"]
        df_gtf = df_gtf[
            df_gtf["len"]
            < int(options.max_exon_len)
        ].copy()
        log.write(
            str(df_gtf.shape[0])
            + "\tlength < max length parameter"
            + os.linesep
        )

        # intersect SUPPA events with this exon list
        df_exons = df_exons[df_exons["ID"].isin(list(df_gtf["name"]))].copy()
        del df_exons["ID"]
        # retain the sorting from SUPPA ioe output file
        # IMPORTANT!
        # suppa in not-deterministic when it comes to the order of compatible transcripts for a given event
        # therefore the output file will be affected
        df_exons.to_csv(options.output, sep="\t", index=False)


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
