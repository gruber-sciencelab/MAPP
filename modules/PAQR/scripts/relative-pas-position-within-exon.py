"""
Calculate relative poly(A) site position within terminal exon
"""

__date__ = "2016-07-07"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument(
    "-v",
    "--verbose",
    dest="verbose",
    action="store_true",
    default=False,
    help="Be loud!",
)

parser.add_argument("--infile", dest="infile", help="TSV table after PAQR infering relative usage")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """
    Main body of the script
    """

    pasPerExon = {}
    with open(options.infile, "r") as infile:
        lines = infile.read().splitlines()
        for line in lines[1:]:
            F = line.rstrip().split("\t")
            if not F[8] in pasPerExon:
                pasPerExon[F[8]] = []
            pasPerExon[F[8]].append(F[3])

    for exon in pasPerExon:
        # sort by pas coordinate
        pas = sorted(pasPerExon[exon])
        # reverse order if the exon is on the negative strand
        pas = pas[::-1] if pas[0].split(":")[2] == "-" else pas
        # keep the original order ~ Maciek
        pas_original_sorting = pasPerExon[exon]

        # get the length of the exon
        name, exon_nr, tot_nr, start, end = exon.split(":")
        start = int(start)
        end = int(end)

        # check if distal poly(A) site is upstream of the end;
        # else: adjust end such that the representative site refers to the end
        chrom, distal_pos, strand = pas[-1].split(":")
        distal_pos = int(distal_pos)
        if strand == "+":
            if end < distal_pos:
                end = distal_pos
        else:
            if start > distal_pos:
                start = distal_pos

        # since exon coordinates are given 1-based, 1 is added to get the distance
        distance = end - start + 1
        for p in pas_original_sorting: # to keep the original order from the input table
            tmp = p.split(":")
            pos = int(tmp[1])
            if strand == "+":
                # get relative position with respect to start coordinate
                sysout("%s\t%.2f\n" % (p, (float(pos - start + 1) / distance * 100.0)))
            else:
                sysout("%s\t%.2f\n" % (p, (float(end - pos + 1) / distance * 100.0)))


if __name__ == "__main__":
    try:
        try:
            options = parser.parse_args()
        except Exception as e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        main(options)
        if options.verbose:
            syserr(
                "### Successfully finished in %i seconds, on %s ###\n"
                % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S"))
            )
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
