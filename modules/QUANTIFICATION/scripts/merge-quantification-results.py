"""
##############################################################################
#
#   Merge salmon quantification results tables of distinct samples
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
        "--quantification-result-directories",
        dest="quant_dirs",
        nargs="+",
        required=True,
        help="Path to the salmon output directory.",
    )
    parser.add_argument(
        "--quantification-genes-filenames",
        dest="gene_quant_file",
        required=True,
        help="File name of the gene-level quantification.",
    )
    parser.add_argument(
        "--quantification-transcripts-filenames",
        dest="transcript_quant_file",
        required=True,
        help="File name of the transcript-level quantification.",
    )
    parser.add_argument(
        "--quantification-merged-genes",
        dest="merged_genes",
        required=True,
        help="Path for the output file with merged genes quantifications.",
    )
    parser.add_argument(
        "--quantification-merged-transcripts",
        dest="merged_transcripts",
        required=True,
        help="Path for the output file with merged transcripts quantifications.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    for infile, outfile in zip(
        [options.gene_quant_file, options.transcript_quant_file],
        [options.merged_genes, options.merged_transcripts],
    ):
        df_list = []
        for res_dir in options.quant_dirs:
            sample_ID = res_dir.split(os.sep)[-1]
            temp = pd.read_csv(os.path.join(res_dir, infile), sep="\t", index_col=0)
            temp = temp[["TPM"]]
            temp.columns = [sample_ID]
            df_list.append(temp)
        merged_df = pd.concat(df_list, axis=1)
        del merged_df.index.name
        merged_df.to_csv(outfile, sep="\t", index_label=False)


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
