"""
##############################################################################
#
#   Create sitecount matrix with raw k-mer counts.
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
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import itertools
import collections
import os


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
        "--sequences",
        dest="fasta",
        required=True,
        help="Path to the input sequences in fasta format.",
    )
    parser.add_argument(
        "--matrix",
        dest="matrix",
        required=True,
        help="Path for the output sitecount matrix file.",
    )
    parser.add_argument(
        "--window-id", dest="window", required=True, help="Encoded window name.",
    )
    parser.add_argument(
        "--kmer-min", dest="kmer_min", required=True, help="Minimal length of kmers."
    )
    parser.add_argument(
        "--kmer-max", dest="kmer_max", required=True, help="Maximal length of kmers."
    )
    return parser


##############################################################################


def count_without_overlaps(seq, sub):
    """
    Generator function to count non-overlapping
    occurances of a sub-string in a string.
    Returns a list of positions at which the substring match starts.
    """
    start = 0
    while True:
        start = seq.find(sub, start)
        if start == -1:
            return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches


def main():
    """Main body of the script."""

    # create a dict of all possible kmers
    all_kmers = {}
    for l in range(int(options.kmer_min), int(options.kmer_max) + 1):
        for i in itertools.product(["A", "C", "G", "T"], repeat=l):
            all_kmers["".join(i)] = 0
    columns_order = sorted(all_kmers.keys())

    # read the fasta input line by line, count kmers, save to the output file:
    with open(options.fasta, "r") as in_file, open(options.matrix, "w") as out_file:

        lines = in_file.read().splitlines()

        # write headers
        out_file.write("coordinates")
        for kmer in columns_order:
            out_file.write("\t" + options.window + "|" + kmer)
        out_file.write(os.linesep)

        # write counts sequence-by-sequence
        for line in range(0, len(lines), 2):
            seq = lines[line + 1]

            # make sure that all the bases are uppercase
            # (input genomic sequence might contain both upper/lower case!)
            seq = seq.upper()

            # In an extremely rare event when a region from reference genome
            # contains N-base we guess the base and have a 25% chance of success.
            # Therefore we guess that all N's are in fact A's:
            seq = seq.replace("N", "A")

            # Here we count kmers:
            #
            # simple count (including k-mer overlaps)
            # (for every sequence the dict gets fully updated)
            # for kmer in columns_order:
            #    all_kmers[kmer] = seq.count(kmer)
            #
            # count without any overlap
            for kmer in columns_order:
                all_kmers[kmer] = len(list(count_without_overlaps(seq, kmer)))

            # save the line with counts
            site_coord = lines[line].split(">")[1]
            out_file.write(site_coord)
            for kmer in columns_order:
                out_file.write("\t" + str(all_kmers[kmer]))
            out_file.write(os.linesep)


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
