"""
##############################################################################
#
#   Merge MotEvo results after scanning a given FASTA file with multiple PWMs.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 18-04-2020
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
import glob
from collections import defaultdict
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
        "--motevo-input-fasta",
        dest="fasta",
        required=True,
        help="Path to the fasta input file for MotEvo runs.",
    )
    parser.add_argument(
        "--motevo-results-dir",
        dest="motevo_dir",
        required=True,
        help="Path to the MotEvo results directory.",
    )
    parser.add_argument(
        "--window-id", dest="window", required=True, help="Encoded window name.",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        required=True,
        help="Path for the output file in TSV format.",
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    # parse MotEvo FASTA input into a list (keep original order)
    # remember to get rid of added MotEvo tag: >>MOTEVO_
    with open(options.fasta) as f:
        fasta_headers = f.read().splitlines()
    fasta_headers = [
        fasta_headers[i].split("_")[-1] for i in range(0, len(fasta_headers), 2)
    ]

    # create a placeholder for the sitecount matrix
    matrix = pd.DataFrame(index=fasta_headers)

    # iterate over all dirs with results per separate motif
    for pwm_dir in glob.glob(options.motevo_dir + "/*"):

        # Open a MotEvo results file.
        # Each record of MotEvo results file consist of 2 lines with the format:
        # [binding position] [strand] [binding posterior] [PWM ID] [binding region]
        # [binding sequence] [binding energy] [FASTA record ID]
        with open(os.path.join(pwm_dir, "posterior_sites")) as f:
            posterior_sites_lines = f.read().splitlines()

        # Iterate over records of MotEvo results.
        # Binding energy should be positive, if it is not then
        # the posterior is probably very low and you should increase the min posterior cutoff.
        # Such records should not be taken into account.
        # Multiple records might correspond to one fasta sequence,
        # if more than one binding site has been found
        # BUT if two sites overlap only the one with the higher posterior is reported by MotEvo!
        binding_dict = defaultdict(list)
        for line in range(0, len(posterior_sites_lines), 2):
            posterior_score = float(posterior_sites_lines[line].split(" ")[2])
            pwm_name = posterior_sites_lines[line].split(" ")[3]
            fasta_ID = posterior_sites_lines[line].split(" ")[4].split("_")[-1]
            energy = float(posterior_sites_lines[line + 1].split(" ")[1])
            if energy > 0.0:
                binding_dict[fasta_ID].append(posterior_score)

        # calculate the summed posterior over all binding positions for a given FASTA record
        summed_posteriors = {k: sum(v) for k, v in binding_dict.items()}
        posteriors_series = pd.Series(summed_posteriors)

        # append the window ID prefix to name a proper column name
        # in the sitecount matrix
        matrix[options.window + "|" + pwm_name] = posteriors_series

    # fill all the missing values with zeros
    matrix = matrix.fillna(0.0)

    # save the sitecount matrix
    matrix.index.rename("coordinates", inplace=True)
    matrix.to_csv(options.outfile, sep="\t")


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
