"""----------------------------------------------------------------------------
This program calculates transcript integrity number (TIN) for each transcript
(or gene) in BED file. TIN is conceptually similar to RIN
(RNA integrity number) but provides transcript level measurement of RNA quality
and is more sensitive to measure low quality RNA samples:

1) TIN score of a transcript is used to measure the RNA integrity
    of the transcript.
2) Median TIN score across all transcripts can be used to measure RNA integrity
    of that "RNA sample".
3) TIN ranges from 0 (the worst) to 100 (the best). TIN = 60 means: 60% of the
    transcript has been covered if the reads coverage were uniform.
4) TIN will be assigned to 0 if the transcript has no coverage or covered reads
    is fewer than cutoff.
----------------------------------------------------------------------------"""
from __future__ import print_function
import math
from multiprocessing import Manager, Pool
from multiprocessing import set_start_method
from optparse import OptionParser
import os
import sys
from time import strftime
import warnings

from bx.intervals import Intersecter, Interval
from qcmodule import BED
import numpy as np
import pysam

warnings.filterwarnings("ignore")

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__ = "2.6.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"
__multithreadingBy__ = "Mihaela Zavolan"


def printlog(mesg):
    """
    print mesg into stderr with time string appending to it.
    """
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    print(mesg, file=sys.stderr)


def uniquefy(seq):
    """
    duplicated members only keep one copy. [1,2,2,3,3,4] => [1,2,3,4].
    """
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def shannon_entropy(arg):
    """
    calculate shannon's H = -sum(P*log(P)). arg is a list of float numbers.
    Note we used natural log here.
    """

    entropy = 0.0

    if not arg:
        return entropy
    # use numpy functions to speed up calculations
    nums = np.array(arg)
    lst_sum = sum(nums)
    fracs = nums / lst_sum
    log_fracs = np.log(fracs)

    entropy = sum(fracs * log_fracs)

    if entropy < 0:
        return -entropy
    return entropy


def build_bitsets(arg_list):
    """
    build intevalTree from list
    """

    ranges = {}
    for element in arg_list:
        chrom = element[0]
        st = element[1]
        end = element[2]
        if chrom not in ranges:
            ranges[chrom] = Intersecter()
        ranges[chrom].add_interval(Interval(st, end))
    return ranges


def union_exons(refbed):
    """
    take the union of all exons defined in refbed file and build bitset
    """

    tmp = BED.ParseBED(refbed)
    all_exons = tmp.getExon()
    unioned_exons = BED.unionBed3(all_exons)
    exon_ranges = build_bitsets(unioned_exons)
    return exon_ranges


def estimate_bg_noise(chrom, tx_st, tx_end, samfile, e_ranges):
    """
    estimate background noise level for a particular transcript
    """

    intron_sig = 0.0  # reads_num * reads_len
    alignedReads = samfile.fetch(chrom, tx_st, tx_end)
    for aligned_read in alignedReads:
        if aligned_read.is_qcfail:
            continue
        if aligned_read.is_unmapped:
            continue
        if aligned_read.is_secondary:
            continue
        read_start = aligned_read.pos
        if read_start < tx_st:
            continue
        if read_start >= tx_end:
            continue
        read_len = aligned_read.qlen
        if len(e_ranges[chrom].find(read_start, read_start + read_len)) > 0:
            continue
        intron_sig += read_len
    return intron_sig


def genomic_positions(refbed, sample_size):
    """
    return genomic positions of each nucleotide in mRNA.
    sample_size: number of nucleotide positions sampled from mRNA.
    """
    if refbed is None:
        print(
            "You must specify a bed file representing gene model\n",
            file=sys.stderr,
        )
        exit(0)

    with open(refbed, "r") as infile:
        for line in infile:
            try:
                if line.startswith(("#", "track", "browser")):
                    continue
                # Parse fields from gene tabls
                fields = line.split()
                chrom = fields[0]
                tx_start = int(fields[1])
                tx_end = int(fields[2])
                geneName = fields[3]
                mRNA_size = sum([int(i) for i in fields[10].strip(",").split(",")])
                exon_starts = [int(x) for x in fields[11].rstrip(",\n").split(",")]
                exon_starts = [x + tx_start for x in exon_starts]
                exon_ends = [int(x) for x in fields[10].rstrip(",\n").split(",")]
                exon_ends = [x + y for (x, y) in zip(exon_starts, exon_ends)]
                intron_size = tx_end - tx_start - mRNA_size
                if intron_size < 0:
                    intron_size = 0
            except Exception:
                print(
                    (
                        "[NOTE:input bed must be 12-column] skipped this "
                        "line: " + line
                    ),
                    file=sys.stderr,
                )
                continue

            chose_bases = [tx_start + 1, tx_end]
            exon_bounds = []
            gene_all_base = []
            if mRNA_size <= sample_size:
                # return all bases of mRNA
                for st, end in zip(exon_starts, exon_ends):
                    # 1-based coordinates on genome, include exon boundaries
                    chose_bases.extend(range(st + 1, end + 1))
                yield (
                    geneName,
                    chrom,
                    tx_start,
                    tx_end,
                    intron_size,
                    chose_bases,
                )
            elif mRNA_size > sample_size:
                step_size = int(mRNA_size / sample_size)
                for st, end in zip(exon_starts, exon_ends):
                    gene_all_base.extend(range(st + 1, end + 1))
                    exon_bounds.append(st + 1)
                    exon_bounds.append(end)
                indx = range(0, len(gene_all_base), step_size)
                chose_bases = [gene_all_base[i] for i in indx]
                yield (
                    geneName,
                    chrom,
                    tx_start,
                    tx_end,
                    intron_size,
                    uniquefy(exon_bounds + chose_bases),
                )


def check_min_reads(samfile, chrom, tx_st, tx_end, cutoff):
    """
    make sure the gene has minimum reads coverage. if cutoff = 10,
    each gene must have 10 *different* reads.
    """

    tmp = False
    read_count = set()
    try:
        alignedReads = samfile.fetch(chrom, tx_st, tx_end)
        for aligned_read in alignedReads:
            if aligned_read.is_qcfail:
                continue
            if aligned_read.is_unmapped:
                continue
            if aligned_read.is_secondary:
                continue
            read_start = aligned_read.pos
            if read_start < tx_st:
                continue
            if read_start >= tx_end:
                continue
            read_count.add(read_start)
            if len(read_count) > cutoff:
                # no need to loop anymore
                tmp = True
                break
        return tmp
    except Exception:
        return False


def genebody_coverage(samfile, chrom, positions, bg_level=0):
    """
    calculate coverage for each nucleotide in *positions*. Sometimes
    len(cvg) < len(positions) because positions where there are no mapped
    reads were ignored.
    """
    cvg = []
    start = positions[0] - 1
    end = positions[-1]
    pos_pnt = -1

    try:
        for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True):
            ref_pos = pileupcolumn.pos + 1
            if ref_pos not in positions:
                continue
            pos_pnt += 1
            # append 0 coverages for positions of interest
            while ref_pos > positions[pos_pnt]:
                pos_pnt += 1
                cvg.append(0.0)
            if pileupcolumn.n == 0:
                cvg.append(0.0)
                continue
            cover_read = 0.0
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del:
                    continue
                if pileupread.alignment.is_qcfail:
                    continue
                if pileupread.alignment.is_secondary:
                    continue
                if pileupread.alignment.is_unmapped:
                    continue
                # if pileupread.alignment.is_duplicate:continue
                cover_read += 1.0
            cvg.append(cover_read)
    except Exception:
        cvg = []

    if bg_level <= 0:
        return cvg
    tmp = []
    for i in cvg:
        subtracted_sig = int(i - bg_level)
        if subtracted_sig > 0:
            tmp.append(subtracted_sig)
        else:
            tmp.append(0)
    return tmp


def tin_score(cvg, length):
    """calculate TIN score"""
    tin = 0
    if len(cvg) == 0 or np.sum(cvg) == 0:
        return tin

    # remove positions with 0 read coverage
    cvg_eff = [float(i) for i in cvg if float(i) > 0]
    entropy = shannon_entropy(cvg_eff)

    tin = 100 * math.exp(entropy) / length

    return tin


# function to run for each gene (data parallelization)
def gf(arg_list):
    (
        sample_TINS_per_transcript,
        gname,
        i_chr,
        i_tx_start,
        i_tx_end,
        intron_size,
        pick_positions,
        sample_name,
        options,
        exon_ranges,
    ) = arg_list

    samfile = pysam.Samfile(sample_name, "rb")

    if gname not in sample_TINS_per_transcript:
        sample_TINS_per_transcript[gname] = []

    noise_level = 0.0

    # check minimum reads coverage
    if (
        check_min_reads(samfile, i_chr, i_tx_start, i_tx_end, options.minimum_coverage)
        is not True
    ):
        sample_TINS_per_transcript[gname] = sample_TINS_per_transcript[gname] + [0.0]
    else:
        # estimate background noise if '-s' was specified
        if options.subtract_bg:
            intron_signals = estimate_bg_noise(
                i_chr, i_tx_start, i_tx_end, samfile, exon_ranges
            )
            if intron_size > 0:
                noise_level = intron_signals / intron_size

        coverage = genebody_coverage(
            samfile, i_chr, sorted(pick_positions), noise_level
        )
        tin1 = tin_score(cvg=coverage, length=len(pick_positions))
        sample_TINS_per_transcript[gname] = sample_TINS_per_transcript[gname] + [tin1]

    samfile.close()
    # log memory usage
    # pid = os.getpid()
    # import psutil
    # py = psutil.Process(pid)
    # memoryUse = py.memory_info()[0]  # memory use in bytes
    # sys.stderr.write('memory use: ' + str(memoryUse) + '\n')


def main():
    set_start_method("spawn")
    usage = "%prog [options]" + "\n" + __doc__ + "\n"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option(
        "-i",
        "--input",
        action="store",
        type="string",
        dest="input_files",
        help=(
            'Input BAM file(s). "-i" takes these input: 1) a single BAM file. '
            '2) "," separated BAM files (no spaces allowed). 3) directory '
            "containing one or more bam files. 4) plain text file containing "
            "the path of one or more bam files (Each row is a BAM file path). "
            "All BAM files should be sorted and indexed using samtools. "
            "[required]"
        ),
    )
    parser.add_option(
        "-r",
        "--refgene",
        action="store",
        type="string",
        dest="ref_gene_model",
        help=(
            "Reference gene model in BED format. Must be strandard 12-column "
            "BED file. [required]"
        ),
    )
    parser.add_option(
        "-c",
        "--minCov",
        action="store",
        type="int",
        dest="minimum_coverage",
        default=10,
        help="Minimum number of read mapped to a transcript. default=%default",
    )
    parser.add_option(
        "-n",
        "--sample-size",
        action="store",
        type="int",
        dest="sample_size",
        default=100,
        help=(
            "Number of equal-spaced nucleotide positions picked from mRNA. "
            "Note: if this number is larger than the length of mRNA (L), it "
            "will be halved until it's smaller than L. default=%default"
        ),
    )
    parser.add_option(
        "--names",
        dest="sample_names",
        action="store",
        type="string",
        help=(
            "sample names, comma separated (no spaces allowed); number must "
            "match the number of provided bam_files"
        ),
    )
    parser.add_option(
        "-s",
        "--subtract-background",
        action="store_true",
        dest="subtract_bg",
        help=(
            "Subtract background noise (estimated from intronic reads). Only "
            "use this option if there are substantial intronic reads."
        ),
    )
    parser.add_option(
        "-p",
        "--processes",
        action="store",
        type="int",
        dest="nrProcesses",
        default=1,
        help=(
            "Number of child processes for the parallelization. Default: 1"
        ),
    )
    (options, args) = parser.parse_args()

    # if '-s' was set
    if options.subtract_bg:
        exon_ranges = union_exons(options.ref_gene_model)
    else:
        exon_ranges = None

    if options.sample_size < 0:
        print("Number of nucleotide can't be negative", file=sys.stderr)
        sys.exit(0)
    elif options.sample_size > 1000:
        print(
            (
                "Warning: '-n' is too large! Please try smaller '-n' value if "
                "program is running slow."
            ),
            file=sys.stderr,
        )

    if not (options.input_files and options.ref_gene_model):
        parser.print_help()
        sys.exit(0)

    if not os.path.exists(options.ref_gene_model):
        print(
            "\n\n" + options.ref_gene_model + " does NOT exist" + "\n",
            file=sys.stderr,
        )
        parser.print_help()
        sys.exit(0)

    printlog("Get BAM file(s) ...")
    bamfiles = options.input_files.split(",")

    if len(bamfiles) <= 0:
        print("No BAM file found, exit.", file=sys.stderr)
        sys.exit(0)
    else:
        print("Total %d BAM file(s):" % len(bamfiles), file=sys.stderr)
        for f in bamfiles:
            print("\t" + f, file=sys.stderr)

    names = options.sample_names.split(",")
    if len(names) != len(bamfiles):
        print(
            "[ERROR] Number of bam files does not match number of names",
            file=sys.stderr,
        )
        sys.exit(2)

    # print header
    sys.stdout.write("transcript")
    for i in names:
        sys.stdout.write("\t%s" % i)
    sys.stdout.write("\n")

    mgr = Manager()
    sample_TINS_per_transcript = mgr.dict()

    genomic_positions_list = genomic_positions(
        refbed=options.ref_gene_model, sample_size=options.sample_size
    )

    for f in bamfiles:
        printlog("Processing " + f)

        conditions = []
        for (
            gname,
            i_chr,
            i_tx_start,
            i_tx_end,
            intron_size,
            pick_positions,
        ) in genomic_positions_list:
            conditions.append(
                [
                    sample_TINS_per_transcript,
                    gname,
                    i_chr,
                    i_tx_start,
                    i_tx_end,
                    intron_size,
                    pick_positions,
                    f,
                    options,
                    exon_ranges,
                ]
            )

        pool = Pool(processes=options.nrProcesses)
        pool.imap_unordered(gf, conditions)

        # clean up processes
        pool.close()
        pool.join()

    for ex in sorted(sample_TINS_per_transcript.keys()):
        vals = [round(x, 10) for x in sample_TINS_per_transcript[ex]]
        print(
            "%s\t%s" % (ex, "\t".join(map(str, vals))),
            file=sys.stdout,
        )


if __name__ == "__main__":
    main()
