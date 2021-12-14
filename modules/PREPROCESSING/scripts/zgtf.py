"""
Downloaded from:
https://github.com/zavolanlab/zgtf/blob/master/zgtf/zgtf.py
(20.06.2020)
Modifications: Maciej Bak
"""

import os
import sys
import operator
import HTSeq

def gtf_to_transcript_exons(gtf, transcript_type):
    """
    Parse gtf and return a dictionary where
    key: trancsript_id
    value: list of exon entries
    """
    gft = HTSeq.GFF_Reader(gtf)

    transcripts = {}

    for gtf_line in gft:
        if gtf_line.type == 'exon':
            try:
                tr_id = gtf_line.attr['transcript_id']
                tr_type = gtf_line.attr['transcript_biotype']
            except:
                sys.stderr.write(f"Problem with: {gtf_line}. Exiting.{os.linesep}")
                sys.exit(1)

            if tr_type != transcript_type:
                continue

            if tr_id not in transcripts:
                transcripts[tr_id] = [gtf_line]
            else:
                transcripts[tr_id].append(gtf_line)

    return transcripts

def transcript_exons_to_bed12(exons_list, transcript_id):
    """
    Convert a list of exon Genomic Intervals from a transcripts (from HTseq) to bed12 line
    """

    blockSizes = []
    blockStarts = []
    sorted_exons = sorted(exons_list, key=operator.attrgetter("iv.start"))
    chrom = sorted_exons[0].iv.chrom
    tr_start = min(sorted_exons[0].iv.start, sorted_exons[0].iv.end, sorted_exons[-1].iv.start, sorted_exons[-1].iv.end)
    tr_end = max(sorted_exons[0].iv.start, sorted_exons[0].iv.end, sorted_exons[-1].iv.start, sorted_exons[-1].iv.end)
    strand = sorted_exons[0].iv.strand
    items = len(sorted_exons)
    
    for exon in sorted_exons:
        blockStarts.append(str(exon.iv.start - tr_start))
        blockSizes.append(str(exon.iv.end - exon.iv.start))
    
    bed12_entry = "\t".join([
        chrom,
        str(tr_start),
        str(tr_end),
        transcript_id,
        "1",
        strand,
        str(tr_start),
        str(tr_end),
        "0",
        str(items),
        ",".join(blockSizes)+",",
        ",".join(blockStarts)+",",
        ])
    
    return bed12_entry