"""Obtain coverage pkl-objects and wiggle files of terminal exon coverage"""

__date__ = "2016-07-11"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import os
import time
import HTSeq
from argparse import ArgumentParser, RawTextHelpFormatter
from collections import defaultdict
import multiprocessing
import itertools
import _pickle as cPickle

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("--bam",
                    dest="bam",
                    help="Name of input bam files to infer the coverage profiles from.")

parser.add_argument("--clusters",
                    dest="polyA",
                    help="file name providing the poly(A) sites in extended BED file format")

parser.add_argument("--min_dist_exStart2prox",
                    dest="min_dist_exStart2prox",
                    type=int, 
                    default=250,
                    help="INT value for the minimum distance between exon start and proximal site. Exon coverage is extended by --upstream_extension if distance is below threshold [default: 250]")

parser.add_argument("--processors",
                    dest="proc",
                    type=int,
                    default=1,
                    help="<INT> maximum number of cores used for multiprocessing [default: 1]")

parser.add_argument("--pickle_out",
                    dest="pkl_out",
                    help="name pickle output file. The name for the wiggle files is inferred from the pickle name (suffix replacement).")

parser.add_argument("--ds_extension", 
                    dest="ds_ext",
                    type=int, 
                    default=200,
                    help="<INT> number of nucleotides appended to the exon end to create coverage also for this region[ default: 200]")

parser.add_argument("--unstranded",
                    dest="unstranded", 
                    action="store_true",
                    default=False,
                    help="Set this option when the data (bam file) is unstranded [default: False]")

# redefine functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def process_single_read(aln,
                        read_regions,
                        region,
                        cvg,
                        splice_var_count,
                        genome_pos_lookup,
                        read_start_of_splice_var,
                        us_exon_start,
                        unstranded):

    strand = region.strand

    splice_candidate = False

    read_weight = 1.0 / aln.optional_field("NH")

    map_types = [ c.type for c in aln.cigar ]
    if "N" in map_types:
        for c_idx in range(1, len( aln.cigar) - 1):
            if aln.cigar[c_idx].type != "N":
                continue
            if ( (strand == "+" and
                  aln.cigar[c_idx - 1].type == "M" and
                  aln.cigar[c_idx + 1].type == "M" and
                  aln.cigar[c_idx - 1].ref_iv.end < region.start and
                  aln.cigar[c_idx].ref_iv.end < region.end and
                  aln.cigar[c_idx + 1].ref_iv.end > region.start) or
                 (strand == "-" and
                  aln.cigar[c_idx + 1].type == "M" and
                  aln.cigar[c_idx - 1].type == "M" and
                  aln.cigar[c_idx + 1].ref_iv.start > region.end and
                  aln.cigar[c_idx].ref_iv.start > region.start and
                  aln.cigar[c_idx - 1].ref_iv.start < region.end) ):
                splice_candidate = True
                break

    if splice_candidate:
        # current read/pair contains a spliced mapper that is used to 
        # get infos about the upstream coverage

        # work only with the mapping parts
        regions_list = []
        for reg in read_regions.steps():
            if "1" in reg[1]:
                regions_list.append( reg[0] )
        border_part_idx = None
        border_read_start = 0
        border_read_end = 0

        if strand == "+":
            # find read part that overlaps with exon start
            for reads_reg_idx in range(len(regions_list)):
                interval = regions_list[ reads_reg_idx ]
                if interval.end >= region.start:
                    # this match overlaps the exon start
                    border_read_start = interval.start
                    border_part_idx = reads_reg_idx
                    break
            # consistency check
            if border_part_idx is None:
                syserr("[ERROR] No part of the current alignment overlaps the exon start border\n")
                syserr("[ERROR] current alignments: %s, %s\n" % (str(aln), str(mate) ) )
                syserr("[ERROR] current matching regions: %s\n" % str(regions_list) )
                return None
            
            # create names for splice variants iteratively
            splice_variant = ''
            # iterate over splice junctions (3' -> 5')
            for exon_idx in range(1,border_part_idx + 1)[::-1]:
                # as coordinates, use really the start/end coordinates here
                # (do not follow: end not included style)
                # for region definition purposes, the end is adjusted later on
                mapped_pos = regions_list[exon_idx].start
                splice_start = regions_list[exon_idx - 1].end - 1
                if splice_variant == '':
                    splice_variant = str(mapped_pos) + ":" + str(splice_start)
                else:
                    splice_variant += ";" + str(mapped_pos) + ":" + str(splice_start)
                
                if splice_variant not in splice_var_count:
                    splice_var_count[ splice_variant] = 1
                    genome_pos_lookup[ splice_variant] = {}
                    us_exon_start[ splice_variant ] = 0
                else:
                    splice_var_count[ splice_variant] += 1

                if exon_idx == 1:
                    # also save the read start of the splice variant
                    if splice_variant not in read_start_of_splice_var:
                        read_start_of_splice_var[ splice_variant ] = []
                    read_start_of_splice_var[ splice_variant ].append( regions_list[0].start )

            # also count all matching parts downstream of the splice border
            for interval in regions_list[border_part_idx:]:
                if((strand == "+" and interval.start > region.end) or
                   (strand == "-" and interval.end < region.start)):
                    continue
                cvg[ interval ] += read_weight

        else:
            # same for negative strand
            # find read part that overlaps with exon end
            # (aka the actual exon start in minus strand direction)
            for reads_reg_idx in range(len(regions_list))[::-1]:
                interval = regions_list[ reads_reg_idx ]
                if interval.start <= region.end:
                    # this match overlaps the exon start
                    border_read_end = interval.end
                    border_part_idx = reads_reg_idx
                    break
            # consistency check
            if border_part_idx is None:
                syserr("[ERROR] No part of the current alignment overlaps the exon end border\n")
                syserr("[ERROR] current alignments: %s, %s\n" % (str(aln), str(mate) ) )
                syserr("[ERROR] current matching regions: %s\n" % str(regions_list) )
                return None

            # crate names for splice variants iteratively
            splice_variant = ''
            # iterate over splice junctions (3' -> 5')
            regions_list_len = len(regions_list)
            for exon_idx in range(border_part_idx, regions_list_len - 1):
                # as coordinates, use really the start/end coordinates here
                # (do not follow: end not included style)
                # for region definition purposes, the end is adjusted later on
                mapped_pos = regions_list[exon_idx].end - 1
                splice_start = regions_list[exon_idx + 1].start
                if splice_variant == '':
                    splice_variant = str(mapped_pos) + ":" + str(splice_start)
                else:
                    splice_variant += ";" + str(mapped_pos) + ":" + str(splice_start)

                if splice_variant not in splice_var_count:
                    splice_var_count[ splice_variant] = 1
                    genome_pos_lookup[ splice_variant] = {}
                    us_exon_start[ splice_variant ] = 0
                else:
                    splice_var_count[ splice_variant] += 1

                if exon_idx == regions_list_len - 2:
                    # also save the read start of the splice variant
                    if splice_variant not in read_start_of_splice_var:
                        read_start_of_splice_var[ splice_variant ] = []
                    read_start_of_splice_var[ splice_variant ].append( regions_list[exon_idx + 1].end - 1 )

            # also count all remaining matching parts
            for interval in regions_list[:border_part_idx + 1]:
                if((strand == "+" and interval.start > region.end) or
                   (strand == "-" and interval.end < region.start)):
                    continue
                cvg[ interval ] += read_weight

    else:
        ##############
        # normal processing of read/pair
                
        # go over read_regions and update cvg vector
        for reg in read_regions.steps():
            # only consider regions that only contain one element
            # '1' was used as element to indicate that this region is valid
            if len( reg[1] ) == 1:
                if(reg[0].start > region.end or
                   reg[0].end < region.start):
                    continue
                cvg[ reg[0] ] += read_weight

    return "finished"

def infer_exon_start(bam, region_string, region_start, region_end, strand, unstranded):

    us_splice_site = 0
    mostUpstream = -1

    # capture the 5' start positions
    # if the reads do not splice out of the region
    # at the 5' end
    five_p_pos = []
    # store the abundance of all splice sites
    # to select the strongest one
    splice_sites = {}

    for aln in bam.fetch(region = region_string):

        # skip alignment if it is not aligned
        if not aln.aligned:
            continue

        # skip alignment if it is on the wrong strand
        if not unstranded:
            # skip alignment if it has the wrong orientation
            if ((aln.pe_which == "first" and aln.iv.strand == strand) or
                (aln.pe_which == "second" and aln.iv.strand != strand)):
                continue
            if (aln.pe_which == 'not_paired_end' and aln.iv.strand != strand):
                    continue

        # make cigar list to be 3' -> 5' oriented
        if strand == "+":
            cigar_list = aln.cigar[::-1]
        else:
            cigar_list = aln.cigar

        # check the 3' most matching part:
        # it's 3' end MUST be upstream of the end of the given region
        # otherwise it can't considered to belong
        # to the selected splice variant
        valid_read = False
        for cig_op in cigar_list:
            if cig_op.type == "M":
                first_match_part = cig_op
                if strand == "+":
                    most_downstream_pos = first_match_part.ref_iv.end - 1
                else:
                    most_downstream_pos = first_match_part.ref_iv.start
                if (most_downstream_pos >= region_start and
                    most_downstream_pos <= region_end):
                    valid_read = True
                break

        if valid_read:
            # check if the read has splice sites
            # and store the upstream pos of the most 3' matching part
            match_regions = []
            for cig_op in aln.cigar:
                if cig_op.type != "M":
                    continue
                match_regions.append( cig_op.ref_iv )
            if len(match_regions) > 1:
                # this read contains a splice-junction
                if strand == "+":
                    splice_pos = match_regions[-1].start
                else:
                    splice_pos = match_regions[0].end - 1
                if splice_pos not in splice_sites:
                    splice_sites[ splice_pos ] = 0
                splice_sites[ splice_pos ] += 1
                    
            else:
                if strand == "+":
                    five_p_pos.append(match_regions[-1].start)
                    if -1 == mostUpstream or mostUpstream > match_regions[-1].start:
                        mostUpstream = match_regions[-1].start
                else:
                    five_p_pos.append(match_regions[0].end)
                    if mostUpstream < match_regions[0].end:
                        mostUpstream = match_regions[0].end

    # get strongest splice variant
    # -> update exon start if necessary
    if len( splice_sites ) > 0:
        top_site = sorted(splice_sites, key=splice_sites.get)[-1]
        top_site_cnt = splice_sites[ top_site ]
        # get the number of reads that extend the top splice site upstream
        if strand == "+":
            longer_noSplice_reads_cnt = sum([1 for x_tmp in five_p_pos if x_tmp < top_site])
        else:
            longer_noSplice_reads_cnt = sum([1 for x_tmp in five_p_pos if x_tmp > top_site])

        if top_site_cnt > longer_noSplice_reads_cnt:
            # splice variant is most supportive
            return (top_site, mostUpstream)
        else:
            return(0, mostUpstream)

    else:
        # if mostUpstream is still set to -1 it means that
        # no valid read was found to extend the exon:
        # hence use previous upstream pos as exon start
        if mostUpstream == -1:
            if strand == "+":
                return( region_start, mostUpstream)
            else:
                return( region_end, mostUpstream)
        return(0, mostUpstream)

def count_same_name_reads(same_read_list,
                          reg_start,
                          reg_end,
                          cvg_tmp,
                          unstranded,
                          strand):

    if len(same_read_list) == 1:
        # count single read alignment
        aln = same_read_list[0]
        read_weight = 1.0 / aln.optional_field("NH")

        if unstranded is True:
            reverse_strand = True if aln.iv.strand != strand else False
        else:
            reverse_strand = True if aln.pe_which == "first" else False
        for cig_op in aln.cigar:
            if cig_op.type != "M":
                continue
            if(cig_op.ref_iv.start > reg_end or
               cig_op.ref_iv.end < reg_start):
                continue
            if reverse_strand:
                cig_op.ref_iv.strand  = "-" if cig_op.ref_iv.strand == "+" else "+"
            cvg_tmp[ cig_op.ref_iv ] += read_weight

    else:
        # read name occurs at least two times
        # find read pairs and process them
        # eventually, process all remaining single end reads
        reiterate = True
        while reiterate:
            mate_name_dict = {}
            mate_pos = None
            for aln_idx in range(len(same_read_list)):

                if aln_idx + 1 == len(same_read_list):
                    # last entry reached
                    # do not iterate over the list again
                    reiterate = False
                            
                curr_aln = same_read_list[aln_idx]
                read_weight = 1.0 / curr_aln.optional_field("NH")

                if not curr_aln.proper_pair:
                    continue

                curr_mate_rank = "first" if curr_aln.pe_which == "second" else "second"
                curr_mate_name = (curr_aln.read.name +
                                  "_" + curr_mate_rank +
                                  ":" +
                                  str(curr_aln.mate_start.start) +
                                  ":" + str( abs(curr_aln.inferred_insert_size) ) )
                if curr_mate_name in mate_name_dict:
                    # pair found
                    # process pair
                    mate_pos = mate_name_dict[ curr_mate_name ]
                    mate = same_read_list [ mate_pos ]

                    # create a genomicArrayOfSets to track the position of the reads
                    # to avoid counting overlapping fragments twice
                    read_regions = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
                    if( (not unstranded and curr_aln.pe_which == "first") or
                        (unstranded and curr_aln.iv.strand != strand)):
                        # current is first of pair (needs strand reversal)
                        # mate is second in pair (correct orientation)
                        for cig_op in curr_aln.cigar:
                            if cig_op.type != "M":
                                continue
                            cig_op.ref_iv.strand = mate.iv.strand
                            read_regions[ cig_op.ref_iv ] += "1"
                        # process mate
                        for cig_op in mate.cigar:
                            if cig_op.type != "M":
                                continue
                            read_regions[ cig_op.ref_iv ] += "1"
                    else:
                        # other case: aln is second in pair, mate first
                        for cig_op in curr_aln.cigar:
                            if cig_op.type != "M":
                                continue
                            read_regions[ cig_op.ref_iv ] += "1"
                        # process mate (and reverse strand)
                        for cig_op in mate.cigar:
                            if cig_op.type != "M":
                                continue
                            cig_op.ref_iv.strand = curr_aln.iv.strand
                            read_regions[ cig_op.ref_iv ] += "1"

                    # go over read_regions and update cvg vector
                    for reg in read_regions.steps():
                        if len( reg[1] ) == 1:
                            if(reg[0].start > reg_end or
                               reg[0].end < reg_start):
                                continue
                            cvg_tmp[ reg[0] ] += read_weight

                    # break the for-loop at the current pos of the list
                    break

                else:
                    # include name into hash
                    curr_name = (curr_aln.read.name +
                                  "_" + curr_aln.pe_which +
                                  ":" +
                                  str(curr_aln.iv.start) +
                                  ":" + str( abs(curr_aln.inferred_insert_size) ) )
                    mate_name_dict[ curr_name ] = aln_idx

            if mate_pos is not None:
                # read pair was processed
                # delete the alns from the list
                same_read_list = [a for i,a in enumerate(same_read_list) if i not in {aln_idx, mate_pos}]

        # iterate over remaining single read alignments
        for a in same_read_list:
            # process a
            read_weight = 1.0 / a.optional_field("NH")
            if unstranded is True:
                reverse_strand = True if a.iv.strand != strand else False
            else:
                reverse_strand = True if a.pe_which == "first" else False
            for cig_op in a.cigar:
                if cig_op.type != "M":
                    continue
                if(cig_op.ref_iv.start > reg_end or
                   cig_op.ref_iv.end < reg_start):
                    continue
                if reverse_strand:
                    cig_op.ref_iv.strand  = "-" if cig_op.ref_iv.strand == "+" else "+"
                cvg_tmp[ cig_op.ref_iv ] += read_weight

    return

def get_coverage(bam,region_string, strand, unstranded):

    # store read pairs and process them together
    curr_lookup = {}
    # store reads with ambiguous mapping
    # (multi mappers for which no clear pair association
    # is possible)
    curr_ambiguous = {}

    cvg_tmp = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' )

    chrom_name, coords = region_string.split(":")
    reg_start, reg_end = [int(x) for x in coords.split("-")]

    aln_list = list(bam.fetch(region = region_string))
    aln_list = sorted(aln_list, key=lambda x: x.read.name)
    same_read_list = []
    for new_aln in aln_list:

        if len(same_read_list) > 0 and new_aln.read.name != same_read_list[0].read.name:
            # process previous read name
            count_same_name_reads(same_read_list,
                                  reg_start,
                                  reg_end,
                                  cvg_tmp,
                                  unstranded,
                                  strand)
            
            # delete previous alignments
            same_read_list = []

        # skip alignment if it is not aligned
        if not new_aln.aligned:
            continue

        if not unstranded:
        # skip alignment if it has the wrong orientation
            if ((new_aln.pe_which == "first" and new_aln.iv.strand == strand) or
                (new_aln.pe_which == "second" and new_aln.iv.strand != strand)):
                continue
            if (new_aln.pe_which == 'not_paired_end' and new_aln.iv.strand != strand):
                    continue

        same_read_list.append(new_aln)

    #----------------------------------------
    # process last read name
    if len(same_read_list) > 0:
        count_same_name_reads(same_read_list,
                              reg_start,
                              reg_end,
                              cvg_tmp,
                              unstranded,
                              strand)
    return cvg_tmp

def process_same_name_reads(same_read_list,
                            region,cvg,
                            splice_var_count,
                            genome_pos_lookup,
                            read_start_of_splice_var,
                            us_exon_start,
                            unstranded,
                            strand):

    if len(same_read_list) == 1:
        # process a solitary alignment
        aln = same_read_list[0]
        if unstranded is True:
            reverse_strand = True if aln.iv.strand != strand else False
        else:
            reverse_strand = True if aln.pe_which == "first" else False
        # create a genomicArrayOfSets for the read
        # this allows similar processing as for read pairs
        read_regions = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
        for cig_op in aln.cigar:
            if cig_op.type != "M":
                continue
            if reverse_strand:
                cig_op.ref_iv.strand = "-" if cig_op.ref_iv.strand == "+" else "+"
            read_regions[ cig_op.ref_iv ] += "1"

        processing_results = process_single_read(aln,
                                                 read_regions,
                                                 region,
                                                 cvg,
                                                 splice_var_count,
                                                 genome_pos_lookup,
                                                 read_start_of_splice_var,
                                                 us_exon_start,
                                                 unstranded)

        # error handling
        if processing_results is None:
            return None

    else:
        # read name occurs at least two times
        # find read pairs and process them
        # eventually, process all remaining single end reads
        reiterate = True
        while reiterate:
            mate_name_dict = {}
            mate_pos = None
            for aln_idx in range(len(same_read_list)):

                if aln_idx + 1 == len(same_read_list):
                    # last entry reached
                    # do not iterate over the list again
                    reiterate = False
                            
                curr_aln = same_read_list[aln_idx]

                if not curr_aln.proper_pair:
                    continue
                
                curr_mate_rank = "first" if curr_aln.pe_which == "second" else "second"
                curr_mate_name = (curr_aln.read.name +
                                  "_" + curr_mate_rank +
                                  ":" +
                                  str(curr_aln.mate_start.start) +
                                  ":" + str( abs(curr_aln.inferred_insert_size) ) )
                if curr_mate_name in mate_name_dict:
                    # pair found
                    # process pair
                    mate_pos = mate_name_dict[ curr_mate_name ]
                    mate = same_read_list [ mate_pos ]

                    # create a genomicArrayOfSets to track the position of the reads
                    # to avoid counting overlapping fragments twice
                    read_regions = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
                    if( (not unstranded and curr_aln.pe_which == "first") or
                        (unstranded and curr_aln.iv.strand != strand)):
                        # current is first of pair (needs strand reversal)
                        # mate is second in pair (correct orientation)
                        for cig_op in curr_aln.cigar:
                            if cig_op.type != "M":
                                continue
                            cig_op.ref_iv.strand = mate.iv.strand
                            read_regions[ cig_op.ref_iv ] += "1"
                        # process mate
                        for cig_op in mate.cigar:
                            if cig_op.type != "M":
                                continue
                            read_regions[ cig_op.ref_iv ] += "1"
                            
                    else:
                        # other case: aln is second in pair, mate first
                        for cig_op in curr_aln.cigar:
                            if cig_op.type != "M":
                                continue
                            read_regions[ cig_op.ref_iv ] += "1"
                        # process mate (and reverse strand)
                        for cig_op in mate.cigar:
                            if cig_op.type != "M":
                                continue
                            cig_op.ref_iv.strand = curr_aln.iv.strand
                            read_regions[ cig_op.ref_iv ] += "1"

                    if unstranded is True:
                        if strand == "+":
                            read_part = curr_aln if curr_aln.iv.start < mate.iv.start else mate
                        else:
                            read_part = curr_aln if curr_aln.iv.start > mate.iv.start else mate
                    else:
                        read_part = curr_aln if curr_aln.pe_which == "second" else mate

                    processing_results = process_single_read(read_part,
                                                             read_regions,
                                                             region,
                                                             cvg,
                                                             splice_var_count,
                                                             genome_pos_lookup,
                                                             read_start_of_splice_var,
                                                             us_exon_start,
                                                             unstranded)

                    # error handling
                    if processing_results is None:
                        return None

                    # break the for-loop at current pos
                    break

                else:
                    # include name into hash
                    curr_name = (curr_aln.read.name +
                                  "_" + curr_aln.pe_which +
                                  ":" +
                                  str(curr_aln.iv.start) +
                                  ":" + str( abs(curr_aln.inferred_insert_size) ) )
                    mate_name_dict[ curr_name ] = aln_idx
                            
            if mate_pos is not None:
                # read pair was processed
                # delete the alns from the list
                same_read_list = [a for i,a in enumerate(same_read_list) if i not in {aln_idx, mate_pos}]

        # iterate over remaining single read alignments
        for a in same_read_list:
            # process a
                            
            if unstranded is True:
                reverse_strand = True if a.iv.strand != strand else False
            else:
                reverse_strand = True if a.pe_which == "first" else False
            # create a genomicArrayOfSets for the read
            # this allows similar processing as for read pairs
            read_regions = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
            for cig_op in a.cigar:
                if cig_op.type != "M":
                    continue
                if reverse_strand:
                    cig_op.ref_iv.strand = "-" if cig_op.ref_iv.strand == "+" else "+"
                read_regions[ cig_op.ref_iv ] += "1"

            processing_results = process_single_read(a,
                                                     read_regions,
                                                     region,
                                                     cvg,
                                                     splice_var_count,
                                                     genome_pos_lookup,
                                                     read_start_of_splice_var,
                                                     us_exon_start,
                                                     unstranded)

            # error handling
            if processing_results is None:
                return None

    return "finished"

def generate_cvg_vector(region, ex_id, to_extend_upstream, bam_file, unstranded):

    # open bam file
    bam=HTSeq.BAM_Reader( bam_file )

    cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' )
    # initialize cvg just to ensure that the chromosome is existent in the saved pkl-file
    # to guarentee processing later when the relative usage is inferred
    cvg[ HTSeq.GenomicPosition(region.chrom, 1, region.strand) ] += 1

    # store first read of read pair
    read_lookup = defaultdict( lambda:"" )

    # store pairs for later processing
    # in the case of upstream extension
    pair_lookup = {}

    # store the mapping of genome positions
    # true position <-> composite exon position without splicing
    genome_pos_lookup = {}
    splice_var_count = {}
    read_start_of_splice_var = {}

    # save for the splice variants, if the upstream exon start was found
    us_exon_start = {}

    region_string = region.chrom + ":" + str(region.start) + "-" + str(region.end)
    # if unstranded is set, the strand info is very important
    # because it is assumed that all reads come from the correct strand
    strand = region.strand

    to_extend_upstream = int(to_extend_upstream)
    if to_extend_upstream > 0:
        # upstream extension of the terminal exon is necessary

        aln_list = list(bam.fetch(region = region_string))
        aln_list = sorted(aln_list, key=lambda x: x.read.name)
        same_read_list = []
        for new_aln in aln_list:

            if len(same_read_list) > 0 and new_aln.read.name != same_read_list[0].read.name:
                # process previous read name
                process_res = process_same_name_reads(same_read_list,
                                                      region,
                                                      cvg,
                                                      splice_var_count,
                                                      genome_pos_lookup,
                                                      read_start_of_splice_var,
                                                      us_exon_start,
                                                      unstranded,
                                                      strand)
                # error handling
                if process_res is None:
                    cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' )
                    cvg[ HTSeq.GenomicPosition(region.chrom, 1, region.strand) ] += 1
                    return (ex_id, cvg, 0)

                # delete previous alignments
                same_read_list = []

            # skip alignment if it is not aligned
            if not new_aln.aligned:
                continue

            if not unstranded:
                # skip alignment if it has the wrong orientation
                if ((new_aln.pe_which == "first" and new_aln.iv.strand == strand) or
                    (new_aln.pe_which == "second" and new_aln.iv.strand != strand)):
                    continue
                if (new_aln.pe_which == 'not_paired_end' and new_aln.iv.strand != strand):
                    continue

            same_read_list.append( new_aln )

        # --------------------------------------------
        # process last read name
        if len(same_read_list) > 0:
            process_res = process_same_name_reads(same_read_list,
                                                  region,
                                                  cvg,
                                                  splice_var_count,
                                                  genome_pos_lookup,
                                                  read_start_of_splice_var,
                                                  us_exon_start,
                                                  unstranded,
                                                  strand)

            # error handling
            if process_res is None:
                # return a minimal coverage array
                cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' )
                cvg[ HTSeq.GenomicPosition(region.chrom, 1, region.strand) ] += 1
                return (ex_id, cvg, 0)

        ########################

        # go through splice-pos association map
        # always only consider the best supported association

        # store splice junctions differently depending on if they relate to:
        # - the original splice border
        # - the first splice junction upstream of the original one
        # - the second     -||-
        # ...
        sorted_splice_junction = []
        top_variants = []
        us_exon_stats = []
        map_to_pos = None
        upstream_pos_tmp = None

        for splice_j_key in splice_var_count:
            splice_nr = len(splice_j_key.split(";"))
            while len(sorted_splice_junction) < splice_nr:
                sorted_splice_junction.append({})
            splice_nr -= 1 # make the index for the splice_nr 0-based
            sorted_splice_junction[ splice_nr ][ splice_j_key ] = splice_var_count[ splice_j_key ]

        # for every splice_part (first, second, ...):
        # only consider the top splice variant
        for splice_part_idx in range(len(sorted_splice_junction)):
            if splice_part_idx == 0:
                # primary splice_site
                curr_splice_var_count = sorted_splice_junction[ splice_part_idx ]
                top_splice_variant = sorted(curr_splice_var_count, key=curr_splice_var_count.get)[-1]
                map_to_pos, most_downstream = [int(tmp) for tmp in top_splice_variant.split(":")]

                # sanity check
                # if the splice-junction is downstream of the proximal poly(A) site:
                # NO extension; return an EMPTY cvg array 
                if( (strand == "+" and map_to_pos > to_extend_upstream) or
                    (strand == "-" and map_to_pos < to_extend_upstream)):
                    syserr("[INFO] splice-junction (%i) is downstream of proximal cleavage site start(%i)\n"
                           % (map_to_pos, to_extend_upstream))
                    syserr("[INFO] Suspect exon %s does not get assigned coverage\n" % ex_id)
                    cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' )
                    cvg[ HTSeq.GenomicPosition(region.chrom, 1, region.strand) ] += 1
                    return (ex_id, cvg, 0)

                # move coordinate into last intron pos
                if strand == "+":
                    map_to_pos -= 1
                else:
                    map_to_pos += 1
                us_exon_stats.append( [most_downstream] )
                top_variants.append(top_splice_variant)

            else:
                # more upstream splice sites
                # only select keys that match previous top splice variant(s)
                prev_top_splice_var = top_variants[-1]
                curr_splice_var_count = {k: sorted_splice_junction[ splice_part_idx ][k] for k in sorted_splice_junction[ splice_part_idx ] if prev_top_splice_var in k}

                # break here if no new splice variant is supported by the so far top splice sites
                if len( curr_splice_var_count ) == 0:
                    break
                top_splice_variant = sorted(curr_splice_var_count, key=curr_splice_var_count.get)[-1]

                prev_exon_start, most_downstream = [int(tmp) for tmp in top_splice_variant.split(";")[-1].split(":")]

                # compare the number of reads that support the top splice site
                # with the number of reads that start in the previous exon and whose start
                # pos is more upstream than the splice site
                #
                # only consider this top splice variant if it is more supported than reads
                # that read through this splice site; in this case: stop looking for more splice sites
                top_splice_variant_count = curr_splice_var_count[ top_splice_variant ]
                if top_variants[-1] in read_start_of_splice_var:
                    if strand == "+":
                        read_starts_count = sum([1 for x_tmp in read_start_of_splice_var[ top_variants[-1] ] if x_tmp < prev_exon_start])
                    else:
                        read_starts_count = sum([1 for x_tmp in read_start_of_splice_var[ top_variants[-1] ] if x_tmp > prev_exon_start])
                else:
                    read_starts_count = 0

                if top_splice_variant_count > read_starts_count:
                    # current splice variant is most supported
                    top_variants.append(top_splice_variant)

                    # update start pos of downstream exon
                    us_exon_stats[-1].append( prev_exon_start )
                    us_exon_stats.append( [most_downstream] )
                else:
                    # more reads start in the current exon than 
                    # reads that support the major splice site
                    # -> do not append additional splice sites

                    # change back the most_downstream site
                    # for the previous exon
                    not_used_tmp, most_downstream = [int(tmp) for tmp in top_variants[-1].split(";")[-1].split(":")]
                    break

        if map_to_pos is None:
            syserr("[INFO] No spliced reads found for region: %s\n" % (region.chrom + ":" + str(region.start) + "-" + str(region.end) ) )
            return (ex_id, cvg, 0)


        ########################
        # extend the most upstream part of the 
        # major splice variant so long till another splice
        # site more upstream is encountered
 
        us_exon_start_pos = 0
        tmp_search_reg_downstream = us_exon_stats[-1][0]
        # also get a rough estimate of the exon region
        # by the starts of the reads -> use the last top_splice variant for this
        if strand == "+":    
            mostUpstream = min( [x_tmp for x_tmp in read_start_of_splice_var[ top_variants[-1] ] ])
        else:
            mostUpstream = max( [x_tmp for x_tmp in read_start_of_splice_var[ top_variants[-1] ] ])

        while 0 == us_exon_start_pos:

            if  strand == "+":
                region_string = region.chrom + ":" + str(mostUpstream) + "-" + str(most_downstream + 1)
                us_value, newUpstream = infer_exon_start(bam,
                                                         region_string,
                                                         mostUpstream,
                                                         tmp_search_reg_downstream,
                                                         strand,
                                                         unstranded)
                tmp_search_reg_downstream = mostUpstream
                mostUpstream = newUpstream
            else:
                region_string = region.chrom + ":" + str(most_downstream) + "-" + str(mostUpstream + 1)
                us_value, newUpstream = infer_exon_start(bam,
                                                          region_string,
                                                          tmp_search_reg_downstream,
                                                          mostUpstream,
                                                          strand,
                                                          unstranded)
                tmp_search_reg_downstream = mostUpstream
                mostUpstream = newUpstream
            us_exon_start_pos = us_value

        # sanity check
        try:
            assert(len(us_exon_stats[-1]) == 1)
        except AssertionError:
            syserr("[ERROR] Most upstream exon already has a start coodinate\n")
            syserr("[ERROR] Most upstream exon: upstream pos: %i, downstream pos: %i\n"
                   % (us_exon_stats[-1][0], us_exon_stats[-1][1]))
            sys.exit(2)

        us_exon_stats[-1].append(us_exon_start_pos)

        ########################
        # create the coverage for the upstream exon(s)

        us_exon_cvgs = []

        for part_idx in range(len(us_exon_stats)):
            ds_pos = us_exon_stats[part_idx][0]
            us_pos = us_exon_stats[part_idx][1]
            if  strand == "+":
                region_string = (region.chrom + ":" + 
                                 str(us_pos) + 
                                 "-" + str(ds_pos + 1) )
            else:
                region_string = (region.chrom + ":" +
                                 str(ds_pos) + "-" +
                                 str(us_pos + 1) )

            us_exon_cvg = get_coverage(bam,region_string, strand, unstranded)
            us_exon_cvgs.append(us_exon_cvg)

        ########################
        # append the upstream exon coverage(s) to the terminal exon
        # use the top splice variant(s) as connection

        length_extension = 0

        for ex_idx in range(len(us_exon_cvgs)):
            ds_pos = us_exon_stats[ex_idx][0]
            us_pos = us_exon_stats[ex_idx][1]
            if strand == "+":
                interval = HTSeq.GenomicInterval( region.chrom, us_pos, ds_pos + 1, "+" )
                us_exon_cvg_list = list(us_exon_cvgs[ex_idx][  interval ].steps())[::-1]
            else:
                interval = HTSeq.GenomicInterval( region.chrom, ds_pos, us_pos + 1, "-" )
                us_exon_cvg_list = list(us_exon_cvgs[ex_idx][  interval ].steps())

            length_extension += interval.length
            for reg in us_exon_cvg_list:
                step_len = reg[0].end - reg[0].start
                if strand == "+":
                    map_interval = HTSeq.GenomicInterval( region.chrom, map_to_pos - step_len + 1, map_to_pos + 1, "+")
                    map_to_pos = map_to_pos - step_len
                else:
                    map_interval = HTSeq.GenomicInterval( region.chrom, map_to_pos, map_to_pos + step_len, "-")
                    map_to_pos = map_to_pos + step_len
                cvg[ map_interval ] += reg[1]
                # test
                #cvg[ reg[0] ] += reg[1]

        return (ex_id, cvg, length_extension)

    # no upstream extension necessary
    else:
        cvg = get_coverage(bam,
                           region_string,
                           strand,
                           unstranded)

        return (ex_id, cvg, 0)

def generate_cvg_vector_star( a_b_c_d ):
    """This is only a wrapper function:
    unpack tuple and invoke the actual function"""
    return generate_cvg_vector( *a_b_c_d )


def main(options):
    """Main logic of the script"""

    # initialize the pool for multimapping
    pool = multiprocessing.Pool( options.proc )

    # check proper suffix of pickle file
    try:
        assert( options.pkl_out.endswith(".pkl") )
    except AssertionError:
        syserr("[ERROR] Output file does not seem to be a correct name for a pickle file (.pkl suffix required)\n")
        sys.exit(2)

    # proprocess option if data is unstranded
    unstranded = options.unstranded

    # process and save poly(A) sites
    polyAsites = defaultdict(list)
    polyA_file = options.polyA
    with open(polyA_file) as f:
            for line in f:
                F = line.strip().split("\t")
                # only consider the first 10 columns for the poly(A) sites
                F = F[0:10]
                polyAsites[ F[8] ].append( F )

    # store exon regions
    exons = []
    # store the exon ids
    exon_ids = []
    # store if the proximal site is close to the exon start
    prox_close2start = []
    adjusted_ex_coords = 0
    for exon in polyAsites:
        # consistency check: does the exon has at least two poly(A) sites
        if len( polyAsites[exon] ) < 2:
            syserr("[ERROR] " + str(exon.split(":")[0]) + "does not have at least two poly(A) sites\n")
            sys.exit(2)

        # parse exon specifications and save the regions
        exon_id, ex_nr, ex_tot, ex_start, ex_end = exon.split(":", 5)
        chrom = polyAsites[ exon ][0][0]
        ex_start = int(ex_start)
        ex_end = int(ex_end)
        # adjust exon coordinates to match
        # BED- (and HTSeq) conventions:
        # 0-based and end excluded
        ex_start -= 1
        strand = polyAsites[ exon ][0][5]
        # extend end of exon if necessary
        distal = polyAsites[ exon ][-1]
        if strand == "+" and int(distal[2]) > ex_end:
            ex_end = int(distal[2])
            adjusted_ex_coords += 1
        elif strand == "-" and int(distal[1]) < ex_start:
            ex_start = int(distal[1])
            adjusted_ex_coords += 1

        # append additional region to the end of the exon, just to generate the coverage profile also there
        if strand == "+":
            ex_end += options.ds_ext
        else:
            ex_start -= options.ds_ext

        # check the most proximal site and it's distance to the exon start
        # store the upstream coordinate of the poly(A) site in case the distance is too small
        if (strand == "+" and
            int(polyAsites[ exon ][0][1]) - ex_start < options.min_dist_exStart2prox):
            prox_close2start.append(polyAsites[ exon ][0][1])
        elif (strand == "-" and
              ex_end - int(polyAsites[ exon ][0][2]) < options.min_dist_exStart2prox):
            prox_close2start.append(polyAsites[ exon ][0][2])
        else:
            # store minus one
            prox_close2start.append(-1)
            
        exons.append(HTSeq.GenomicInterval( chrom, ex_start, ex_end, strand ) )
        exon_ids.append(exon)

    syserr("[INFO] Number of exon for which the exon end was extended because of the poly(A) site coordinates: %i\n" % (adjusted_ex_coords))

    ###
    ### Finished first part
    ###
    ### now iterate over bam files and over regions

    syserr("[INFO] bam file to process:\n")
    syserr("[INFO] %s\n\n" % (options.bam))

    out_file_basename = options.pkl_out.replace(".pkl", "")
    exon_extensions = {}

    # store coverage
    cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' )

    # process regions in prallel
    cvg_results_tuple_list = pool.map(generate_cvg_vector_star, 
                                      zip(exons,
                                          exon_ids,
                                          prox_close2start,
                                          itertools.repeat(options.bam), 
                                          itertools.repeat(unstranded)
                                      )
    )
    # Merge all coverage parts into one coverage file
    for res_tuple in cvg_results_tuple_list:
        if res_tuple is None:
            continue
        ex_id = res_tuple[0]
        cvg_part = res_tuple[1]
        extension = res_tuple[2]
        exon_extensions[ ex_id ] = extension
        for c in cvg_part.steps():
            cvg[ c[0] ] += c[1]

    # # test
    # for ex_idx in range(len(exons)):
    #     cvg_part = generate_cvg_vector_star((exons[ex_idx], prox_close2start[ex_idx], options.bam, unstranded))
    #     for c in cvg_part.steps():
    #         cvg[ c[0] ] += c[1]

    # write coverage files
    with open(options.pkl_out, 'wb') as output:
        cPickle.dump(cvg, output, -1)

    # write wiggle files
    cvg.write_bedgraph_file( out_file_basename + ".plus.wig", "+" )
    cvg.write_bedgraph_file( out_file_basename + ".minus.wig", "-" )

    # write exon-specific extensions
    with open( out_file_basename + ".extensions.tsv", "w") as ext_out:
        for ex in exon_extensions:
            ext_out.write("%s\t%i\n" % (ex, exon_extensions[ex]))


if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception as e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)
        pkl_path = "/".join( options.pkl_out.split("/")[:-1])
        if pkl_path != "":
            # check of the output directory for the wiggle files
            # already exists; otherwise, make the directory
            if not os.path.isdir(pkl_path):
                # create the directory
                # (do not handle a potential OSError)
                os.makedirs(pkl_path)
                syserr("[INFO] Output directory for pickle objects created\n")
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
