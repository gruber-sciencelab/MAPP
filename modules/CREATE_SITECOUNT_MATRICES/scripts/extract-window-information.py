"""
##############################################################################
#
#   A wrapper for tasks related to extraction of window-wise information:
#   Given a list of sites in bed format, genome in fasta and an encoded
#   relative coordinates of windows -> extract absolute window coordinates 
#   and their sequences.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 09-04-2020
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
import pybedtools


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
        "--window-id",
        dest="window",
        required=True,
        help="Encoded relative coordinates of a genomic window.",
    )
    parser.add_argument(
        "--bed",
        dest="bed",
        required=True,
        help="Path to a bed-formatted file with the coordinates of sites of interest.",
    )
    parser.add_argument(
        "--genome",
        dest="genome",
        required=True,
        help="Path to a genome file in fasta format.",
    )
    parser.add_argument(
        "--outdir", dest="outdir", required=True, help="Path for the output directory."
    )
    return parser


##############################################################################


def main():
    """Main body of the script."""

    # test that the input BED-file indeed contains coordinates for single sites
    with open(options.bed) as bed_file:
        bed_lines = bed_file.read().splitlines()
    for line in bed_lines:
        line_parsed = line.split("\t")
        assert int(line_parsed[2]) - int(line_parsed[1]) == 1

    outdir = options.outdir

    # parse the relative coordinates from the encoded window name
    coords = [i[1:].split("to") for i in options.window.split(".")]
    coords = [coords[0][0], coords[0][1], coords[1][0], coords[1][1]]

    # Construct bash commands to extract region coordinates into bed
    if int(coords[0]) and not int(coords[3]):  # whole window upstream the site
        assert not int(coords[2])
        beg = coords[0]
        end = coords[1]
        command = (
            "cat "
            + options.bed
            + ' | awk -F "\\t" \
        \'{{ if ($6 == "+") print $1"\\t"$2-'
            + beg
            + '"\\t"$2-'
            + end
            + '"\\t"$4"\\t"$5"\\t"$6; else if ($6 == "-") print $1"\\t"$3+'
            + end
            + '"\\t"$3+'
            + beg
            + '"\\t"$4"\\t"$5"\\t"$6}}\' \
        1> '
            + os.path.join(outdir, "coordinates.bed")
        )

    elif int(coords[3]) and not int(coords[0]):  # whole window downstream the site
        assert not int(coords[1])
        beg = coords[2]
        end = coords[3]
        command = (
            "cat "
            + options.bed
            + ' | awk -F "\\t" \
        \'{{ if ($6 == "+") print $1"\\t"$2+'
            + beg
            + '"\\t"$2+'
            + end
            + '"\\t"$4"\\t"$5"\\t"$6; else if ($6 == "-") print $1"\\t"$3-'
            + end
            + '"\\t"$3-'
            + beg
            + '"\\t"$4"\\t"$5"\\t"$6}}\' \
        1> '
            + os.path.join(outdir, "coordinates.bed")
        )

    else:  # window goes through the site
        assert not int(coords[1]) and not int(coords[2])
        beg = coords[0]
        end = coords[3]
        command = (
            "cat "
            + options.bed
            + ' | awk -F "\\t" \
        \'{{ if ($6 == "+") print $1"\\t"$2-'
            + beg
            + '"\\t"$2+'
            + end
            + '"\\t"$4"\\t"$5"\\t"$6; else if ($6 == "-") print $1"\\t"$3-'
            + end
            + '"\\t"$3+'
            + beg
            + '"\\t"$4"\\t"$5"\\t"$6}}\' \
        1> '
            + os.path.join(outdir, "coordinates.bed")
        )

    # call bash cat->awk command to extract absolute coordinates of the window
    os.system(command)

    # prevent errors due to filesystem latency (comp. cluster issues)
    time.sleep(60)

    # Extract sequences under windows coordinates
    BedTool = pybedtools.BedTool(os.path.join(outdir, "coordinates.bed"))
    BedTool = BedTool.sequence(fi=options.genome, s=True)
    with open(os.path.join(outdir, "sequences.fasta"), "w") as outfile:
        fasta_lines = open(BedTool.seqfn).read().splitlines()
        assert len(fasta_lines) / 2 == len(bed_lines)  # sanity check
        for i in range(0, len(fasta_lines), 2):
            bed_coordinate = bed_lines[int(i / 2)].split("\t")[3]
            outfile.write(">" + bed_coordinate + os.linesep)
            outfile.write(fasta_lines[i + 1] + os.linesep)


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
