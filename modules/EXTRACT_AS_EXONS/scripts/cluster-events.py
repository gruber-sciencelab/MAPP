"""
##############################################################################
#
#   Cluster the set of alternatively spliced exons
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
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
from collections import defaultdict, namedtuple
import pandas as pd
import scipy.cluster


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
        "--events",
        dest="events",
        required=True,
        help="Path to the file with skipped exon events.",
    )
    parser.add_argument(
        "--outfile",
        dest="outfile",
        required=True,
        help="Path for the file with clustered events.",
    )
    parser.add_argument(
        "--max-distance",
        dest="maxdist",
        required=True,
        help="Maximal distance for the objects with one cluster.",
    )
    return parser


class SkippedExonEvent:
    """
    Implements Skipped Exon (SE) event from SUPPA
    """

    def __init__(self, ID, genes, supporting_transcripts, compatible_transcripts):
        """takes: string,list,list"""
        self.genes = genes
        self.name = ID
        self.supporting = supporting_transcripts
        self.compatible = compatible_transcripts
        self.chrom = ID.split(":")[0]
        self.exon_range = ID.split(":")[1]
        self.beg = int(self.exon_range.split("-")[0])
        self.end = int(self.exon_range.split("-")[1])
        self.strand = ID.split(":")[-1]


def get_event_ID(ev):
    """
    extracts unique ID of a skipped exon event
    """
    chrom = ev.split(";")[1][3:].split(":")[0]
    return chrom + ":" + ev.split("-")[1].replace(":", "-") + ":" + ev.split(":")[-1]


def min_coverage_distance(e1_beg, e1_end, e2_beg, e2_end):
    """
    returns [0,1]-normalized overlap distance between two exons
    """
    Range = namedtuple("Range", ["start", "end", "length"])
    r1 = Range(start=e1_beg, end=e1_end, length=e1_end - e1_beg + 1)
    r2 = Range(start=e2_beg, end=e2_end, length=e2_end - e2_beg + 1)
    latest_start = max(r1.start, r2.start)
    earliest_end = min(r1.end, r2.end)
    overlap = earliest_end - latest_start + 1  # inclusive ranges!
    if overlap <= 0:
        return 1
    return 1 - min(overlap / r1.length, overlap / r2.length)


##############################################################################


def main():
    """Main body of the script."""

    SUPPA_e = pd.read_csv(options.events, sep="\t")
    SUPPA_e["ID"] = SUPPA_e.apply(
        lambda row: get_event_ID(row["event_id"]), axis=1
    )  # add unique IDs

    # parse the data into a dictionary
    unified_events = {}
    for i, row in SUPPA_e.iterrows():
        if row.ID in unified_events.keys():
            unified_events[row.ID][0] = unified_events[row.ID][0] + "," + row.gene_id
            unified_events[row.ID][1] = (
                unified_events[row.ID][1] + "," + row.alternative_transcripts
            )
            unified_events[row.ID][2] = (
                unified_events[row.ID][2] + "," + row.total_transcripts
            )
        else:
            unified_events[row.ID] = [
                row.gene_id,
                row.alternative_transcripts,
                row.total_transcripts,
            ]

    # construct SkippedExonEvent objects
    event_list = []
    for ev in unified_events.keys():
        genes = list(set(unified_events[ev][0].split(",")))
        genes = ",".join(sorted(genes))  # this will be an unique wellordered key
        included = list(set(unified_events[ev][1].split(",")))
        compatible = list(set(unified_events[ev][2].split(",")))
        event_list.append(SkippedExonEvent(ev, genes, included, compatible))

    # group the events by genes
    gene_dict = defaultdict(list)
    for ev in event_list:
        gene_dict[ev.genes].append(ev)

    # create a dataframe with SE events
    collapsed_events_columns = [
        "genes",
        "chr",
        "start",
        "end",
        "strand",
        "included",
        "total",
    ]
    collapsed_events = pd.DataFrame(columns=collapsed_events_columns)
    for g in gene_dict.keys():
        for ev in gene_dict[g]:
            # remove suppa non-deterministic order of transcripts here
            supporting = ",".join(sorted(ev.supporting))
            compatible = ",".join(sorted(ev.compatible))
            row = [g, ev.chrom, ev.beg, ev.end, ev.strand, supporting, compatible]
            collapsed_events.loc[len(collapsed_events)] = row

    ### CLUSTERING ###
    cluster_index = 0
    with open(options.outfile, "w") as f:
        events_by_chromosomes = collapsed_events.groupby(
            "chr"
        )  # first group events by chromosomes
        for name, chromosome in events_by_chromosomes:

            plus_strand = chromosome[chromosome["strand"] == "+"].reset_index()
            if len(plus_strand) == 0:  # extreme case, watch out for contig chromosomes
                pass
            elif (
                len(plus_strand) == 1
            ):  # extreme case,  watch out for contig chromosomes
                ID = (
                    plus_strand.at[0, "chr"]
                    + ":"
                    + str(plus_strand.at[0, "start"])
                    + "-"
                    + str(plus_strand.at[0, "end"])
                    + ":"
                    + str(plus_strand.at[0, "strand"])
                )
                f.write(
                    ID
                    + "\t"
                    + plus_strand.at[0, "genes"]
                    + "\t"
                    + plus_strand.at[0, "chr"]
                    + "\t"
                    + str(plus_strand.at[0, "start"])
                    + "\t"
                    + str(plus_strand.at[0, "end"])
                    + "\t"
                    + plus_strand.at[0, "strand"]
                    + "\t"
                )
                f.write(
                    plus_strand.at[0, "included"]
                    + "\t"
                    + plus_strand.at[0, "total"]
                    + "\t"
                    + str(cluster_index)
                    + os.linesep
                )
                cluster_index += 1
            else:
                # cluster those on + strand:
                # create a condensed distance matrix as in 'scipy'
                condensed_distances = []
                for i in range(0, len(plus_strand)):
                    for j in range(i + 1, len(plus_strand)):
                        val = min_coverage_distance(
                            plus_strand.at[i, "start"],
                            plus_strand.at[i, "end"],
                            plus_strand.at[j, "start"],
                            plus_strand.at[j, "end"],
                        )
                        condensed_distances.append(val)
                # and cluster based on the distance matrix:
                linkage = scipy.cluster.hierarchy.linkage(
                    y=condensed_distances, method="complete",
                )
                clusters = scipy.cluster.hierarchy.fcluster(
                    Z=linkage, t=float(options.maxdist), criterion="distance"
                )
                # iterate over clusters
                for label in list(set(clusters)):
                    indices = [i for i, l in enumerate(list(clusters)) if l == label]
                    # iterate over all the exons in this clusters and save all the events
                    for i, row in plus_strand.iloc[indices].iterrows():
                        ID = (
                            str(row.chr)
                            + ":"
                            + str(row.start)
                            + "-"
                            + str(row.end)
                            + ":"
                            + str(row.strand)
                        )
                        f.write(
                            ID
                            + "\t"
                            + row.genes
                            + "\t"
                            + str(row.chr)
                            + "\t"
                            + str(row.start)
                            + "\t"
                            + str(row.end)
                            + "\t"
                            + row.strand
                            + "\t"
                        )
                        f.write(
                            row.included
                            + "\t"
                            + row.total
                            + "\t"
                            + str(cluster_index)
                            + os.linesep
                        )
                    cluster_index += 1

            minus_strand = chromosome[chromosome["strand"] == "-"].reset_index()
            if len(minus_strand) == 0:  # extreme case, watch out for contig chromosomes
                pass
            elif (
                len(minus_strand) == 1
            ):  # extreme case, watch out for contig chromosomes
                ID = (
                    minus_strand.at[0, "chr"]
                    + ":"
                    + str(minus_strand.at[0, "start"])
                    + "-"
                    + str(minus_strand.at[0, "end"])
                    + ":"
                    + str(minus_strand.at[0, "strand"])
                )
                f.write(
                    ID
                    + "\t"
                    + minus_strand.at[0, "genes"]
                    + "\t"
                    + minus_strand.at[0, "chr"]
                    + "\t"
                    + str(minus_strand.at[0, "start"])
                    + "\t"
                    + str(minus_strand.at[0, "end"])
                    + "\t"
                    + minus_strand.at[0, "strand"]
                    + "\t"
                )
                f.write(
                    minus_strand.at[0, "included"]
                    + "\t"
                    + minus_strand.at[0, "total"]
                    + "\t"
                    + str(cluster_index)
                    + os.linesep
                )
                cluster_index += 1
            else:
                # cluster those on - strand
                # create a condensed distance matrix as in 'scipy'
                condensed_distances = []
                for i in range(0, len(minus_strand)):
                    for j in range(i + 1, len(minus_strand)):
                        val = min_coverage_distance(
                            minus_strand.at[i, "start"],
                            minus_strand.at[i, "end"],
                            minus_strand.at[j, "start"],
                            minus_strand.at[j, "end"],
                        )
                        condensed_distances.append(val)
                # and cluster based on the distance matrix:
                linkage = scipy.cluster.hierarchy.linkage(
                    y=condensed_distances, method="complete",
                )
                clusters = scipy.cluster.hierarchy.fcluster(
                    Z=linkage, t=float(options.maxdist), criterion="distance"
                )
                # iterate over clusters
                for label in list(set(clusters)):
                    indices = [i for i, l in enumerate(list(clusters)) if l == label]
                    # iterate over all the exons in this clusters and save all the events
                    for i, row in minus_strand.iloc[indices].iterrows():
                        ID = (
                            str(row.chr)
                            + ":"
                            + str(row.start)
                            + "-"
                            + str(row.end)
                            + ":"
                            + str(row.strand)
                        )
                        f.write(
                            ID
                            + "\t"
                            + row.genes
                            + "\t"
                            + str(row.chr)
                            + "\t"
                            + str(row.start)
                            + "\t"
                            + str(row.end)
                            + "\t"
                            + row.strand
                            + "\t"
                        )
                        f.write(
                            row.included
                            + "\t"
                            + row.total
                            + "\t"
                            + str(cluster_index)
                            + os.linesep
                        )
                    cluster_index += 1


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
