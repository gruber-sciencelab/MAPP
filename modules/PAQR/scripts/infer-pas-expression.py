# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:30:36 2016

Define poly(A) sites and calculate relative usages parallelized

@author: schmiral
"""
__date__ = "Wed Oct 05 13:30:36 2016"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"


# imports
import sys
import os
import time
import HTSeq
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import _pickle as cPickle
import datetime
import multiprocessing
import bisect

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
                    
parser.add_argument("--clusters",
                    dest="polyA",
                    help=("file name providing the poly(A) " +
                        "sites in extended BED file format") )

parser.add_argument("--coverages",
                    dest="coverages",
                    nargs="*",
                    help=("name of the pickle files that " +
                        "store the HTSeq objects with coverages") )

parser.add_argument("--conditions",
                    dest="conditions",
                    nargs="*",
                    help=("values that match the coverages " +
                        "parameters and indicate the condition each " +
                        "coverage object belongs to") )

parser.add_argument("--design_file",
                    dest="design_file",
                    required=True,
                    help=("Path to the tsv table with quality-filtered samples)") )
                    
parser.add_argument("--read_length", 
                    dest="read_length",
                    type=int,
                    default=100,
                    help=("integer value for the expected read length " +
                        "(used for window to fit a regression line to the " +
                        "upstream region of poly(A) sites) [default: 100]") )

parser.add_argument("--min_coverage_region",
                    dest="min_coverage_region",
                    type=int, 
                    default=100,
                    help=("integer value that defines the minimum region " +
                        "that needs to be available to calculate a mean " +
                        "coverage and the mean squared error [default: 100]") )

parser.add_argument("--min_mean_coverage",
                    dest="min_mean_coverage",
                    type=float, 
                    default=20.0,
                    help=("float value that defines the minimum mean " +
                        "coverage required for an exon to be considered " +
                        "[default: 20.0]") )

parser.add_argument("--min_cluster_distance",
                    dest="min_cluster_distance",
                    type=int, 
                    default=200,
                    help=("integer value that defines the distance till " +
                        "which clusters are merged and analyzed separately " +
                        "[default: 200]") )
                    
parser.add_argument("--mse_ratio_threshold",
                    dest="mse_ratio_threshold",
                    type=float,
                    default=0.7,
                    help=("the ratio of the mse after including a new break " +
                        "point devided by the mse before must be below this " +
                        "threshold [default: 0.7]") )
                    
parser.add_argument("--best_break_point_upstream_extension",
                    dest="best_break_point_upstream_extension",
                    type=int,
                    default=200,
                    help=("nucleotide extension added to the cluster's " +
                        "upstream end during the check for global break " +
                        "points [default: 200]") )
                    
parser.add_argument("--ds_reg_for_no_coverage",
                    dest="ds_reg_for_no_coverage",
                    type=int, 
                    default=200,
                    help=("integer value that defines the region that is " +
                        "downstream of the exon and necessary to search for " +
                        "zero coverage downstream of the distal site  [default: 200]") )

parser.add_argument("--max_downstream_coverage",
                    dest="max_downstream_coverage",
                    default=5,
                    type=float,
                    help=("<FLOAT> maximum percentage of the coverage at the exon " +
                          "start is allowed for the mean coverage in the downstream " +
                          "region of a valid distal site [default: 5]") )

parser.add_argument("--processors",
                    dest="proc",
                    type=int,
                    default=1,
                    help="<INT> Number of processors used [default: 1]")

parser.add_argument("-d",
                    "--debug",
                    dest="debug",
                    action="store_true",
                    default=False,
                    help="Print additional diagnostic comments for single exons")

parser.add_argument("--ex_extension_files",
                    dest="ex_extension_files",
                    nargs="*",
                    help=("name of the files that " +
                        "store the upstream extensions per exon for the corresponding sample") )

parser.add_argument("--expressions_out",
                    dest="expressions_out",
                    help="name of file to write the raw expression values per site per sample to")

parser.add_argument("--distal_sites",
                    dest="distal_sites",
                    help="file name to which the exons are written that have only distal site usage\n")

parser.add_argument("--names",
                    dest="sample_names",
                    nargs="*",
                    help="names of the samples")

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

np.seterr(all='raise')

def time_now():#return time
    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")            
    
def find_distal_site(site,
                     ex_start,
                     strand, 
                     region, 
                     merged_exons,
                     max_downstream_coverage,
                     read_length, 
                     ds_reg_extend,
                     debug):
    
    """Check if the passed poly(A)site is applicable as distal site
    
    Use the combined coverage of all samples here"""

    #----------------------------------------

    # define the exon start reference region
    # get the mean coverage of this region
    # (ex_start saves the strand specific exon start)
    if strand == "+":
        ex_start_idx = region.index(ex_start)

        if ex_start + read_length > site[1].start:
            # consistency check:
            # is exon start more upstream than site start
            if site[1].start < ex_start:
                syserr("[INFO] Distal site search stopped whith the current site %s that overlaps the exon start\n"
                       % site[0] )
                return(None, -1)

            # exon is very short: use the start of the currently investigated
            # poly(A) site as boundary
            ex_start_2_idx = region.index(site[1].start)
        else:
            ex_start_2_idx = region.index(ex_start + read_length)
    else:
        if ex_start - read_length < site[1].end:
            # consistency check:
            # is exon start more upstream than site start (meaning the strand specific starts)
            if site[1].end >= ex_start:
                syserr("[INFO] Distal site search stopped whith the current site %s that overlaps the exon start\n"
                       % site[0] )
                return(None, -1)
            
            # exon is very short: use the start of the currently investigated
            # poly(A) site as boundary
            ex_start_idx = region.index( site[1].end)
        else:
            ex_start_idx = region.index(ex_start - read_length)
            
        ex_start_2_idx = region.index(ex_start - 1) + 1

    if ex_start_2_idx - ex_start_idx <= 0:
        syserr("[INFO] Distal site search stopped whith the current site %s that overlaps the exon start\n"
               % site[0])
        return(None, -1)
        
    ex_start_cov = np.mean(merged_exons[ ex_start_idx:ex_start_2_idx ])

    if ex_start_cov == 0:
        if debug:
            syserr( ("[INFO] No coverage at exon start of site %s\n" +
                     "[INFO] Found during search for distal site. Exon is skipped\n")
                    % site[0] )
        return (None, -1)

    #----------------------------------------

    # define the upstream reference region
    # get the mean coverage of this region

    if strand == "+":
        if site[1].start - (2 * read_length) < ex_start:
            # exon is very short
            # use the exon start as ref-region start
            ex_us_idx = region.index(ex_start)
        else:
            ex_us_idx = region.index( site[1].start - (2 * read_length) )
        ex_us_2_idx = region.index( site[1].start )
    else:
        ex_us_idx = region.index( site[1].end )
        if site[1].end + (2 * read_length) >= ex_start:
            # exon is very short
            # use the exon start as ref-region start
            ex_us_2_idx = region.index(ex_start - 1) + 1
        else:
            ex_us_2_idx = region.index( site[1].end + (2 * read_length) )
    ex_us_cov = np.mean(merged_exons[ ex_us_idx:ex_us_2_idx ])

    #----------------------------------------

    # define the downstream reference region
    # get the mean coverage of this region

    if strand == "+":
        ex_ds_idx = region.index(site[1].end )
        # ex_ds_2_idx = region.index( site[1].end + (2 * read_length) )
        ex_ds_2_idx = region.index( site[1].end + (ds_reg_extend) - 1 ) + 1
        ex_ds_cov = np.mean( merged_exons[ ex_ds_idx:ex_ds_2_idx] )
    else:
        # ex_ds_idx = region.index( site[1].start - (2 * read_length) )
        ex_ds_idx = region.index( site[1].start - (ds_reg_extend) )
        ex_ds_2_idx = region.index( site[1].start )
        ex_ds_cov = np.mean(merged_exons[ ex_ds_idx:ex_ds_2_idx ])

    #----------------------------------------

    # a valid distal site has a higher us_mean_cvg than ds_mean_cvg
    # and the percentage of the ds_mean_cvg with respect to
    # the start_mean_cvg is below the threshold
    # (threshold defines the maximum percentage the ds_region is allowed to have)

    if( ex_us_cov > ex_ds_cov and
        ex_ds_cov / float(ex_start_cov) * 100 <= max_downstream_coverage):
        # distal site found
        return(1, 1)
    else:
        # do not use this site as distal site
        return (None, 0)

def assess_minimum_cvg(exon_structure, bp_coverages_dict, region, min_mean_cvg):
    '''Iterate over the dict with coverages and assess
    for every single sample if it has sufficient mean covearge
    at the current exon
    '''

    return_dict = {}
    
    start = exon_structure[2]
    end = exon_structure[3]

    # infer the indices for the current exon boundaries
    start_idx = region.index( start )
    try:
        end_idx = region.index( end )
    except ValueError:
        # for negative strand sites:
        # exon end may not be accessible as position
        end_idx = region.index( end - 1) + 1

    for cond in bp_coverages_dict:
        return_dict[ cond ] = []
        for cvg in bp_coverages_dict[cond]:
            if np.mean( cvg[start_idx:end_idx] ) < min_mean_cvg:
                return_dict[ cond ].append(0)
            else:
                return_dict[ cond ].append(1)

    return return_dict

def find_best_break_point_in_full_exon(bp_cvg_dict,
                                       valid_mean_cvg_dict,
                                       region,
                                       strand,
                                       ex_start,
                                       ex_end,
                                       min_cvg_reg,
                                       mse_ratio_threshold,
                                       upstream_extension):
    '''this function is only used when only a distal site
    was found and the upstream coverage is screened for
    unknown sites'''

    # store for which samples a valid break point was inferred
    valid_break_point_dict = {}

    # prepare the arrays for the search of the best separation point
    curr_ex_start_dict = {}
    curr_ex_end_dict = {}
    # infer sample-specific exon boundaries
    # (this is necessary due to the individual upstream extension)
    for cond in bp_cvg_dict:
        curr_ex_start_dict[cond]=[]
        curr_ex_end_dict[cond]=[]
        for cvg_idx in range(len(bp_cvg_dict[cond])):
            # extend the exon start by the upstream extension
            curr_add_us_ext = upstream_extension[cond][cvg_idx]
            if strand == "+":
                # use the extended exon as start pos
                curr_ex_start_dict[cond].append( ex_start - curr_add_us_ext )
                curr_ex_end_dict[cond].append( ex_end )
            else:
                # minus strand
                # extend the exon at the other end
                    curr_ex_start_dict[cond].append( ex_start )
                    curr_ex_end_dict[cond].append( ex_end + curr_add_us_ext )

    for cond in bp_cvg_dict:
        valid_break_point_dict[cond] = []
        for cvg_idx in range(len(bp_cvg_dict[cond])):

            cvg = bp_cvg_dict[cond][cvg_idx]

            # if the current sample has no valid mean cvg:
            # store False in valid_break_point_dict
            if valid_mean_cvg_dict[cond][cvg_idx] == 0:
                valid_break_point_dict[cond].append(False)
                continue

            # initialize the best mse with stupidly high value
            curr_best_ratio = 10

            # set of indices for the region
            try:
                upstream_start_idx = region.index( curr_ex_start_dict[cond][cvg_idx] )
            except ValueError:
                syserr("[ERROR] Upstream start pos not in list while scanning for overall break point\n")
                syserr("[ERROR] Upstream start: %i\n" 
                       %( curr_ex_start_dict[cond][cvg_idx] ) )
                raise ValueError("[ERROR] Interrupt due to missing position\n")

            try:
                downstream_end_idx = region.index( curr_ex_end_dict[cond][cvg_idx] - 1 ) + 1
            except ValueError:
                syserr("[ERROR] Downstream end pos not in list\n")
                syserr("[ERROR] Downstream end: %i\n" 
                       %( curr_ex_end_dict[cond][cvg_idx] ) )
                raise ValueError("[ERROR] Interrupt due to missing position\n")

            # global mse
            cvg_global = cvg[upstream_start_idx:downstream_end_idx]
            global_mse = np.mean( (cvg_global - np.mean(cvg_global))**2 )

            best_pos_idx = None
            for pos_idx in range(upstream_start_idx + min_cvg_reg + 1, downstream_end_idx - min_cvg_reg - 1):
                # get mse for the current two parts
                cvg_right = cvg[pos_idx + 1: downstream_end_idx]
                cvg_left = cvg[upstream_start_idx:pos_idx]
                cov_diff = cvg_left - np.mean(cvg_left)
                cov_diff = np.append( cov_diff, (cvg_right - np.mean(cvg_right) ) )
                curr_mse = np.mean( cov_diff**2 )
                mse_ratio = float(curr_mse) / global_mse
                if mse_ratio < curr_best_ratio:
                    curr_best_ratio = mse_ratio

            if curr_best_ratio < mse_ratio_threshold:
                valid_break_point_dict[cond].append(True)
            else:
                valid_break_point_dict[cond].append(False)

    return valid_break_point_dict
    
def check_mse_improvement(bp_coverages_dict,
                          valid_mean_cvg_dict,
                          region, 
                          site, 
                          strand, 
                          segment_start,
                          segment_end,
                          min_coverage_region, 
                          best_break_point_upstream_extension, 
                          mse_ratio_threshold, 
                          upstream_extension, 
                          prev_mse,
                          extended_prev_mse=None):
    """check if current poly(A) site improves the mean squared error"""
    
    us_mse_dict = {}
    ds_mse_dict = {}
    mse_ratio_dict = {}

    segment_starts = {}
    segment_ends = {}

    # for both strands: iterate over region from 3' to 5'
    fit_region_to_iterate = []
    fit_region_start = site[1].start
    fit_region_end = site[1].end

    if strand == "+" :
        fit_region_start -= best_break_point_upstream_extension
        fit_region_to_iterate = range( fit_region_start, fit_region_end )[::-1]
    else:
        fit_region_end += best_break_point_upstream_extension
        fit_region_to_iterate = range( fit_region_start, fit_region_end )

    for cond in bp_coverages_dict:
        us_mse_dict[cond] = []
        ds_mse_dict[cond] = []
        mse_ratio_dict[cond] = []
        segment_starts[cond] = []
        segment_ends[cond] = []
        
        cvg_idx = -1
        for cvg in bp_coverages_dict[cond]:
            cvg_idx += 1
            curr_best_ratio = 10

            # store whether the exon was already extended umstream
            # (artificial exon upstream extension only possible once per site)
            already_extended = False

            # if the current sample has no valid mean cvg:
            # append 10 as curr_best_ratio and continue
            if valid_mean_cvg_dict[cond][cvg_idx] == 0:
                mse_ratio_dict[cond].append(curr_best_ratio)
                # append a placeholder value also as segment start/end (strand dependent)
                if strand == "+":
                    segment_starts[cond].append(0)
                else:
                    segment_ends[cond].append(0)
                continue

            # set up up- and downstream region
            # (depending on the strand, some values are "None" still afterwards)
            if strand == "+":
                up_start = segment_start
                down_end = segment_end
                # store the upstream starts in case these are extended in the next for loop
                segment_starts[cond].append(up_start)
            else:
                down_start = segment_start
                up_end = segment_end
                # store the upstream end (aka strand specific starts) in case these are extended in the next for loop
                segment_ends[cond].append(up_end)
                
            prev_reg_mse = prev_mse[cond][cvg_idx]
            if prev_reg_mse == 0:
                # skip this coverage because for this region it is already 
                # completely flat
                syserr( ("[INFO] Skip cond %s cvg %i. Found previous reg " + 
                       "to have zero mse for site %s!\n")
                       % (cond, cvg_idx, str(site[0])) )
                # no reasonable position was found: still the starting value of 10
                # is appended as curr_best_ratio
                mse_ratio_dict[cond].append(curr_best_ratio)
                continue

            for pos in fit_region_to_iterate:
                
                # leave out the current pos (
                # pos does not belong to the upstream nor to the downstream region)
                if strand == "+":
                    up_end = pos
                    down_start = pos + 1
                else:
                    up_start = pos + 1
                    down_end = pos

                if(up_end - 
                    up_start < min_coverage_region):
                    # upstream region is too small
                    # if currently the proximal site is considered
                    # the additional upstream extension of coverage is exploited
                    # current site is the most proximal == extended_prev_mse != None
                    if extended_prev_mse is not None and not already_extended:
                        # extended upstream region is used, hence, 
                        # the extended previous mse should be used as well
                        prev_reg_mse = extended_prev_mse[cond][cvg_idx]

                        if strand == "+":
                            up_start -= upstream_extension[cond][cvg_idx]
                            segment_starts[cond][cvg_idx] = up_start
                        else:
                            up_end += upstream_extension[cond][cvg_idx]
                            segment_ends[cond][cvg_idx]= up_end
                        already_extended = True

                        if(up_end - 
                            up_start < min_coverage_region):
                            # syserr(("[INFO] No upstream region for %s at pos " +
                            #         "%i and potential more upstream positions " +
                            #         "even though upstream extension of the exon. " +
                            #         "Skip this and more upstream positions.\n") 
                            #         % (site[0], pos))
                            break
                    else:
                        # syserr(("[INFO] No upstream region for internal pA %s at pos %i. " +
                        #         "Skip this and more upstream positions.\n")
                        #         % (site[0], pos))
                        break
                if(down_end - 
                    down_start <= min_coverage_region):
                    syserr("[INFO] No downstream region for %s at pos %i\n" 
                            % (site[0], pos))
                    break

                ### Finished definition of indexes

                ### set up upstream region
                try:
                    up_start_idx = region.index(up_start)
                except ValueError:
                    syserr("[ERROR] upstream start pos not in list\n")
                    syserr("[ERROR] Upstream start: %i, upstream stop: %i, diff: %i, min_diff: %i\n" 
                           %( up_start, up_end, up_end - up_start, min_coverage_region) )
                    raise ValueError("[ERROR] Interrupt due to missing position\n")

                try:
                    up_end_idx = region.index(up_end)
                except ValueError:
                    # for negative strand sites:
                    # exon end not may not be accessible as position
                    up_end_idx = region.index(up_end - 1) + 1

                cvg_upstream = cvg[up_start_idx:up_end_idx]
                mean_upstream = np.mean( cvg_upstream )

                ### set up downstream region
                down_start_idx = region.index(down_start)
                down_end_idx = region.index(down_end)

                cvg_downstream = cvg[down_start_idx:down_end_idx]
                mean_downstream = np.mean( cvg_downstream )

                # continue if the downstream mean is higher than the upstream mean
                if mean_downstream > mean_upstream:
                    continue

                # get a mean squared error of up- and downstream region
                cov_diff = cvg_downstream - mean_downstream
                cov_diff = np.append( cov_diff, cvg_upstream - mean_upstream )
                curr_mse = np.mean( cov_diff**2)
                
                mse_ratio = float(curr_mse) / prev_reg_mse
                if mse_ratio < curr_best_ratio:
                    curr_best_ratio = mse_ratio

            # if no reasonable position was found: still the starting value of 10
            # is appended as curr_best_ratio
            mse_ratio_dict[cond].append(curr_best_ratio)
            
    # check if site has valid mse improvement in any of the conditions
    # for this the median mse ratio must meet the threshold criterion
    valid_conds = []
    for cond in mse_ratio_dict:
        # only consider samples with valid mean_cvg
        mse_ratio_subset = np.array([mse_ratio_dict[cond][i] for i in range(len( mse_ratio_dict[cond] ) ) if valid_mean_cvg_dict[cond][i] == 1])
        if len(mse_ratio_subset) == 0:
            continue
        
        # different possibilities a condition is considered as "valid":
        # two samples or less: both must be below the threshold
        # three to five samples: the median must be below the threshold
        # six or more samples: at least three samples must be below th threshold
        if len(mse_ratio_subset) <= 2:
            if sum( mse_ratio_subset < mse_ratio_threshold ) == len(mse_ratio_subset):
                valid_conds.append(cond)
        elif len(mse_ratio_subset) <= 5:
            if len(mse_ratio_subset) % 2 == 1 and np.median( mse_ratio_subset ) < mse_ratio_threshold:
                valid_conds.append(cond)
            elif len(mse_ratio_subset) % 2 == 0 and np.median( np.append(mse_ratio_subset,0) ) < mse_ratio_threshold:
                # the number of samples is even, therefore a zero is appended temporarily
                # this ensures that the median corresponds to half of the samples
                valid_conds.append(cond)
        else:
            if sum( mse_ratio_subset < mse_ratio_threshold) >= 3:
                valid_conds.append(cond)
            
    if len(valid_conds) > 0:
        
        # potential site found
        # complete infos for mse dicts
        # define the new up_mse and ds_mse values based on 
        # the poly(A) site coordinates
        for cond in bp_coverages_dict:
            cvg_idx = -1
            for cvg in bp_coverages_dict[cond]:

                cvg_idx += 1

                if valid_mean_cvg_dict[cond][cvg_idx] == 0:
                    us_mse_dict[cond].append(100)
                    ds_mse_dict[cond].append(100)
                    continue

                # calculate the missing up- and downstream mses for the coverages
                # exclude the region of the poly(A) site cluster when calculating
                # the mses
                if strand == "+":
                    up_start_idx = region.index(segment_starts[cond][cvg_idx])
                    up_end_idx = region.index(site[1].start)
                    down_start_idx = region.index(site[1].end)
                    down_end_idx = region.index(segment_end)
                else:
                    up_start_idx = region.index(site[1].end)
                    try:
                        up_end_idx = region.index(segment_ends[cond][cvg_idx])
                    except ValueError:
                        # for negative strand sites:
                        # exon end may not be accessible as position
                        up_end_idx = region.index(segment_ends[cond][cvg_idx] - 1) + 1
                    down_start_idx = region.index(segment_start)
                    down_end_idx = region.index(site[1].start)

                cvg_upstream = cvg[up_start_idx:up_end_idx]
                mean_upstream = np.mean(cvg_upstream)
                us_mse = np.mean((cvg_upstream - mean_upstream)**2)

                cvg_downstream = cvg[down_start_idx:down_end_idx]
                mean_downstream = np.mean( cvg_downstream )
                ds_mse = np.mean((cvg_downstream - mean_downstream)**2)

                us_mse_dict[cond].append(us_mse)
                ds_mse_dict[cond].append(ds_mse)

    return (mse_ratio_dict, us_mse_dict, ds_mse_dict, valid_conds)

def check_global_best_breaking_point(bp_cvg,
                                     used_sites_positions,
                                     start_idx,
                                     end_idx,
                                     min_coverage_region):

    """Find best breaking point irrespective of poly(A) sites.
    Function to process single sites."""

    cvg_global = bp_cvg[start_idx:end_idx]
    global_mse = np.mean( (cvg_global - np.mean(cvg_global))**2 )
    smallest_mse = global_mse
    best_pos_idx = None
    for pos_idx in range(start_idx + min_coverage_region, end_idx - min_coverage_region):
        # get mse for the current two parts
        cvg_right = bp_cvg[pos_idx + 1: end_idx]
        cvg_left = bp_cvg[start_idx:pos_idx]
        cov_diff = cvg_left - np.mean(cvg_left)
        cov_diff = np.append( cov_diff, (cvg_right - np.mean(cvg_right) ) )
        curr_mse = np.mean( cov_diff**2 )
        if curr_mse < smallest_mse:
            smallest_mse = curr_mse
            best_pos_idx = pos_idx

    if best_pos_idx is not None:
                
        if best_pos_idx not in used_sites_positions:
            # break point is outside of regions that are poly(A)-site supported
            return(0, best_pos_idx)

        else:
            # valid break point found
            return( best_pos_idx, best_pos_idx )

    # no break point was found; return "None"
    return( best_pos_idx, best_pos_idx )
    
def check_global_best_breaking_points(bp_cvg,
                                      cond,
                                      sample_idx,
                                      region,
                                      used_sites_positions,
                                      still_undetected,
                                      start_idx,
                                      end_idx,
                                      min_coverage_region):
    """Find best breaking points irrespective of poly(A) sites.
    Wrapper for the entire exon."""
    
    boundary_indices = [start_idx, end_idx]
    while any(x is not None for x in still_undetected):
        check_curr_segment = False
        for segment_idx in range(1,len(boundary_indices)):
            for idx in range(boundary_indices[segment_idx -1], boundary_indices[segment_idx]):
                if idx in used_sites_positions and still_undetected[used_sites_positions[idx]] is not None:
                    check_curr_segment = True
                    break
            if check_curr_segment:
                break
        
        # use the current borders to find the best breaking point
        (pos_idx_indicator, pos_idx) = check_global_best_breaking_point(bp_cvg,
                                                                        used_sites_positions,
                                                                        boundary_indices[segment_idx -1],
                                                                        boundary_indices[segment_idx],
                                                                        min_coverage_region)
        if pos_idx_indicator == 0:
            # syserr("[ERROR] Break point outside of poly(A) site regions\n")
            return (boundary_indices, pos_idx)
        elif pos_idx_indicator is None:
            # syserr("[ERROR] No valid break point found")
            return (boundary_indices, pos_idx)

        else:
            bisect.insort(boundary_indices, pos_idx)
            still_undetected[used_sites_positions[pos_idx]] = None
    
    return "Finished"
    
def find_sites_with_usage(sites,
                          segment_start,
                          segment_end,
                          segment_mse,
                          original_extended_segment_mse,
                          bp_coverages_dict,
                          valid_mean_cvg_dict,
                          region,
                          strand,
                          min_coverage_region,
                          best_break_point_upstream_extension, 
                          mse_ratio_threshold, 
                          add_us_extension):
    """iterate over the poly(A) sites
    and find the best poly(A) site
    """

    # keep the original value of extended_segment_mse untouched
    # to know  whether it needs to be updated below in the code
    extended_segment_mse = original_extended_segment_mse

    top_site = None
    top_site_idx = None
    top_site_prot_support = 0
    top_mse_ratio = 1
    top_us_mse_dict = {}
    top_ds_mse_dict = {}
    top_site_valid_conditions = []
    top_site_supporting_cvgs = {}
    top_site_is_from_cluster = None
    
    for site_idx in range(len(sites)):
        site = sites[site_idx]

        if isinstance(site[0], list):
            # set the extended_segment_mse to None if the most proximal
            # site has enough distance to the exon start
            # (extended_segment_mse might be None anyways already 
            #  if the current segment is a subsegment of the exon that
            #  is not most proximally located)
            if( (strand == "+" and (site[0][1].start
                                    - segment_start
                                    - best_break_point_upstream_extension) >= min_coverage_region) or
                (strand == "-" and (segment_end
                                    - site[0][1].end
                                    - best_break_point_upstream_extension) >= min_coverage_region)):
                extended_segment_mse = None

            # process a set of closely spaced clusters
            best_candidate_ratio = 1
            best_candidate_prot_support = 0
            best_candidate_up_mse_dict = {}
            best_candidate_ds_mse_dict = {}
            best_candidate_site = None
            best_candidate_valid_conditions = []
            best_candidate_supporting_cvgs = {}
            # iterate over set of closely spaced sites (prox -> dist):
            subsite_idx = -1
            for subsite in site:
                subsite_idx += 1
                (mse_ratios, 
                 us_mse_dict,
                 ds_mse_dict,
                 valid_conds_list) = check_mse_improvement(bp_coverages_dict,
                                                           valid_mean_cvg_dict,
                                                           region, 
                                                           subsite, 
                                                           strand, 
                                                           segment_start,
                                                           segment_end,
                                                           min_coverage_region,  
                                                           best_break_point_upstream_extension, 
                                                           mse_ratio_threshold, 
                                                           add_us_extension,
                                                           segment_mse,
                                                           extended_segment_mse)
                 
                # Use the MINIMAL median-ratio of this site (only from samples below the cutoff)
                # check all conditions separately
                # when ratio is the same: take the site with higher protocol support
                # when two sites share the same top ratio AND the protocol support:
                # take the more proximal site
                curr_best_median_ratio = 2
                curr_supporting_cvgs = {}
                for cond in valid_conds_list:
                    valid_ratios = [ i for i in mse_ratios[cond] if i < mse_ratio_threshold ]
                    curr_supporting_cvgs[cond] = [ i < mse_ratio_threshold for i in mse_ratios[cond]]
                    if np.median( valid_ratios ) < curr_best_median_ratio:
                        curr_best_median_ratio = np.median( valid_ratios )
                if curr_best_median_ratio < best_candidate_ratio:
                    best_candidate_ratio = curr_best_median_ratio
                    best_candidate_site = subsite
                    best_candidate_prot_support = subsite[2]
                    best_candidate_up_mse_dict = us_mse_dict
                    best_candidate_ds_mse_dict = ds_mse_dict
                    best_candidate_valid_conditions = valid_conds_list
                    best_candidate_supporting_cvgs = curr_supporting_cvgs
                elif( curr_best_median_ratio == best_candidate_ratio and
                      subsite[2] > best_candidate_prot_support):
                    best_candidate_ratio = curr_best_median_ratio
                    best_candidate_site = subsite
                    best_candidate_prot_support = subsite[2]
                    best_candidate_up_mse_dict = us_mse_dict
                    best_candidate_ds_mse_dict = ds_mse_dict
                    best_candidate_valid_conditions = valid_conds_list
                    best_candidate_supporting_cvgs = curr_supporting_cvgs
            if best_candidate_ratio < 1:
                # for this set of sites, a potentially used sites was found
                if best_candidate_ratio < top_mse_ratio:
                    top_site = best_candidate_site
                    top_site_idx = site_idx
                    top_mse_ratio = best_candidate_ratio
                    top_site_prot_support = best_candidate_prot_support
                    top_us_mse_dict = best_candidate_up_mse_dict
                    top_ds_mse_dict = best_candidate_ds_mse_dict
                    top_site_valid_conditions = best_candidate_valid_conditions
                    top_site_supporting_cvgs = best_candidate_supporting_cvgs
                    top_site_is_from_cluster = subsite_idx
                    
        else:
            # site is solitary

            # only use the extended_segment_mse if necessary
            if( (strand == "+" and (sites[site_idx][1].start
                                    - segment_start
                                    - best_break_point_upstream_extension) >= min_coverage_region) or
                (strand == "-" and (segment_end
                                    - site[1].end
                                    - best_break_point_upstream_extension) >= min_coverage_region)):
                extended_segment_mse = None

            site = sites[site_idx]
            (mse_ratios, 
             us_mse_dict,
             ds_mse_dict,
             valid_conds_list) = check_mse_improvement(bp_coverages_dict,
                                                       valid_mean_cvg_dict,
                                                       region,
                                                       site,
                                                       strand,
                                                       segment_start,
                                                       segment_end,
                                                       min_coverage_region,
                                                       best_break_point_upstream_extension,
                                                       mse_ratio_threshold,
                                                       add_us_extension,
                                                       segment_mse,
                                                       extended_segment_mse)
            
            curr_best_median_ratio = 2
            curr_supporting_cvgs = {}
            for cond in valid_conds_list:
                valid_ratios = [ i for i in mse_ratios[cond] if i < mse_ratio_threshold ]
                curr_supporting_cvgs[cond] = [ i < mse_ratio_threshold for i in mse_ratios[cond]]
                if np.median( valid_ratios ) < curr_best_median_ratio:
                    curr_best_median_ratio = np.median( valid_ratios )
            if curr_best_median_ratio < top_mse_ratio:
                top_site = site
                top_site_idx = site_idx
                top_mse_ratio = curr_best_median_ratio
                top_site_prot_support = site[2]
                top_us_mse_dict = us_mse_dict
                top_ds_mse_dict = ds_mse_dict
                top_site_valid_conditions = valid_conds_list
                top_site_supporting_cvgs = curr_supporting_cvgs
                top_site_is_from_cluster = None
                
    if top_site is not None:
        # site with usage found

        # update information for this site:
        # valid conditions
        top_site.append(top_site_valid_conditions)
        # samples that indicate usage of this site
        top_site.append(top_site_supporting_cvgs)
        
        upstream_sites = sites[:top_site_idx]
        if len(upstream_sites) > 0:
            # determine the characteristics 
            # of the region upstream of the top site
            if strand == "+":
                up_seg_start = segment_start
                # check if top site belonged to a set of closely spaced clusters
                if top_site_is_from_cluster is None:
                    up_seg_end = sites[top_site_idx][1].start
                else:
                    up_seg_end = sites[top_site_idx][top_site_is_from_cluster][1].start
            else:
                if top_site_is_from_cluster is None:
                    up_seg_start = sites[top_site_idx][1].end
                else:
                    up_seg_start = sites[top_site_idx][top_site_is_from_cluster][1].end
                up_seg_end = segment_end

            if original_extended_segment_mse is not None:
                # update the extended_segment_mse
                extended_segment_mse = original_extended_segment_mse
                for cond in bp_coverages_dict:
                    extended_segment_mse[cond] = []
                    for cvg_idx in range(len(bp_coverages_dict[cond])):
                        cvg = bp_coverages_dict[cond][cvg_idx]
                        if strand == "+":
                            tmp_start_idx =  region.index(segment_start - add_us_extension[cond][cvg_idx])
                            if top_site_is_from_cluster is None:
                                tmp_end_idx = region.index(sites[top_site_idx][1].start)
                            else:
                                tmp_end_idx = region.index(sites[top_site_idx][top_site_is_from_cluster][1].start)
                        else:
                            if top_site_is_from_cluster is None:
                                tmp_start_idx = region.index(sites[top_site_idx][1].end)
                            else:
                                tmp_start_idx = region.index(sites[top_site_idx][top_site_is_from_cluster][1].end)
                            try:
                                tmp_end_idx = region.index(segment_end + add_us_extension[cond][cvg_idx])
                            except ValueError:
                                # the end of the extended exon might not be accessible as position anymore
                                tmp_end_idx = region.index(segment_end + add_us_extension[cond][cvg_idx] - 1) + 1

                        extended_y = cvg[tmp_start_idx:tmp_end_idx]
                        ext_seg_mean = np.mean( extended_y )
                        ext_seg_diff = extended_y - ext_seg_mean
                        ext_seg_mse = np.mean( ext_seg_diff**2)
                        extended_segment_mse[cond].append( ext_seg_mse)
                    
            else:
                # no extended_segment_mse necessary
                extended_segment_mse = None       
                
        if top_site_idx < len(sites) - 2:
            # downstream region also still contains potential sites
            downstream_sites = sites[top_site_idx + 1:]
            if strand == "+":
                if top_site_is_from_cluster is None:
                    ds_seg_start = sites[top_site_idx][1].end
                else:
                    ds_seg_start = sites[top_site_idx][top_site_is_from_cluster][1].end
                ds_seg_end = segment_end
            else:
                ds_seg_start = segment_start
                if top_site_is_from_cluster is None:
                    ds_seg_end = sites[top_site_idx][1].start
                else:
                    ds_seg_end = sites[top_site_idx][top_site_is_from_cluster][1].start

            if len(upstream_sites) > 0:
                # up- and downstream remain sites that need to be checked
                return( find_sites_with_usage(upstream_sites,
                                              up_seg_start,
                                              up_seg_end,
                                              top_us_mse_dict,
                                              extended_segment_mse,
                                              bp_coverages_dict,
                                              valid_mean_cvg_dict,
                                              region,
                                              strand,
                                              min_coverage_region,
                                              best_break_point_upstream_extension, 
                                              mse_ratio_threshold, 
                                              add_us_extension) + 
                        [top_site] + 
                        find_sites_with_usage(downstream_sites,
                                              ds_seg_start,
                                              ds_seg_end,
                                              top_ds_mse_dict,
                                              None, # do not provide an extended upstream mse
                                              bp_coverages_dict,
                                              valid_mean_cvg_dict,
                                              region,
                                              strand,
                                              min_coverage_region,
                                              best_break_point_upstream_extension, 
                                              mse_ratio_threshold, 
                                              add_us_extension))
            else:
                # only downstream remain sites that need to be checked
                return( [top_site] + 
                        find_sites_with_usage(downstream_sites,
                                              ds_seg_start,
                                              ds_seg_end,
                                              top_ds_mse_dict,
                                              None,
                                              bp_coverages_dict,
                                              valid_mean_cvg_dict,
                                              region,
                                              strand,
                                              min_coverage_region,
                                              best_break_point_upstream_extension, 
                                              mse_ratio_threshold, 
                                              add_us_extension))
            
        elif len(upstream_sites) > 0:
            # only upstream remain sites that need to be checked
            return( find_sites_with_usage(upstream_sites,
                                          up_seg_start,
                                          up_seg_end,
                                          top_us_mse_dict,
                                          extended_segment_mse,
                                          bp_coverages_dict,
                                          valid_mean_cvg_dict,
                                          region,
                                          strand,
                                          min_coverage_region,
                                          best_break_point_upstream_extension, 
                                          mse_ratio_threshold, 
                                          add_us_extension) + 
                        [top_site] )
                    
        else:
            # return the current site as only possible site with usage
            return [top_site]
            
    else:
        # no site with usage found anymore
        return []
    
def get_pA_set(exon_id, 
               exon_structure, 
               region, 
               bp_coverages_dict, 
               merged_exons, 
               polyA_sites,
               exon_ext_dict,
               options):
    """Define a set of used poly(A) sites"""

    strand = exon_structure[1]
    # save the strand specific exon start (used during distal site determination)
    ex_start = exon_structure[2] if strand == "+" else exon_structure[3]
    # store the poly(A) sites that are inferred to be used in the current exon
    used_sites = []
    # get the exon upstream extensions
    add_us_ext = exon_ext_dict
    
    ######
    # first part
    ######
    
    ### define the distal site

    distal_site = None
    best_fit_pos = None
    # iterate over the sites from distal to proximal
    for site_idx in range(len(polyA_sites))[::-1]:
            if isinstance(polyA_sites[site_idx][0], list):
                # current "site" is a set of closely spaced sites
                # that are processed together
                site_found = False
                # iterate over the set of sites, from distal to proximal
                for subsite_idx in range(len(polyA_sites[site_idx]))[::-1]:
                    site = polyA_sites[site_idx][subsite_idx]
                    (best_fit_pos, best_fit) = find_distal_site(site,
                                                                ex_start,
                                                                strand,
                                                                region,
                                                                merged_exons, 
                                                                options.max_downstream_coverage,
                                                                options.read_length,
                                                                options.ds_reg_for_no_coverage,
                                                                options.debug)
                    if best_fit == -1:
                        # skip exon because it has no exon start coverage
                        # or the currently processed site already overlaps
                        # with the exon start
                        site_found = True
                        break
                        
                    if best_fit_pos is not None:
                        # a valid regression line was fitted to 
                        # the current upstream region
                        distal_site = site
                        site_found = True
                        break
                if site_found:
                    # break also outer loop that runs over all poly(A)sites
                    break
                    
            else:
                # the current site is a solitary site
                site = polyA_sites[site_idx]
                (best_fit_pos, best_fit) = find_distal_site(site,
                                                            ex_start,
                                                            strand,
                                                            region,
                                                            merged_exons,
                                                            options.max_downstream_coverage,
                                                            options.read_length,
                                                            options.ds_reg_for_no_coverage,
                                                            options.debug)
                if best_fit == -1:
                    # skip exon because it has no exon start coverage
                    # or the currently processed site already overlaps
                    # with the exon start
                    break
                    
                if best_fit_pos is not None:
                    distal_site = site
                    # a valid regression line was fitted to 
                    # the current upstream region
                    break

    if distal_site is None:
        if options.debug:
            syserr("[INFO] No distal site found for %s. Exon skipped\n" 
                   % exon_id)
        return ([], {})
        
    # limit the set of sites to the ones upstream of the distal site
    # (if the distal sites comes from a set of closely spaced sites
    # the entire set is marked as distal and excluded from further analysis)
    polyA_sites = polyA_sites[0:site_idx]
        
    # adjust the exon end according to the distal site end
    if strand == "+":
        exon_structure[3] = distal_site[1].end
    else:
        exon_structure[2] = distal_site[1].start

    ######
    # second part
    ######

    # check if the mean coverage per sample is above the given cutoff
    # return a dict with the conditions as keys and a list as value
    # that holds 0 or 1 at idx-positions that correspond to the coverages
    valid_mean_cvg_dict = assess_minimum_cvg(exon_structure, bp_coverages_dict, region, options.min_mean_coverage)

    if sum([i for lst in valid_mean_cvg_dict.values() for i in lst]) == 0:
        if options.debug:
            syserr("[INFO] No sample with sufficient coverage for %s. Exon skipped\n" 
                   % exon_id)
        return ([],valid_mean_cvg_dict)

    ######
    # third part
    ######

    # handle cases in which only the distal sites remained
    if len(polyA_sites) == 0:
        # check exon region for global best break point
        curr_ex_start = exon_structure[2]
        curr_ex_end = exon_structure[3]
        valid_break_point_dict = find_best_break_point_in_full_exon(bp_coverages_dict,
                                                                    valid_mean_cvg_dict,
                                                                    region,
                                                                    strand,
                                                                    curr_ex_start,
                                                                    curr_ex_end,
                                                                    options.min_coverage_region,
                                                                    options.mse_ratio_threshold,
                                                                    add_us_ext)
        
        
        # return distal site together with the valid_break_point_dict!
        # NOT the valid_mean_cvg_dict!
        return (polyA_sites, valid_break_point_dict)

    ######
    # fourth part
    ######
            
    # find sites with usage
    
    # define the mse of the total region for all samples
    mse_total_region = {}
    mse_total_ext_region = {}
    for cond in bp_coverages_dict:
        mse_total_region[cond] = []
        mse_total_ext_region[cond] = []
        cvg_idx = -1
        for cvg in bp_coverages_dict[cond]:
            cvg_idx += 1

            # do not consider coverages with little cvg
            if valid_mean_cvg_dict[cond][cvg_idx] == 0:
                mse_total_region[cond].append(0)
                mse_total_ext_region[cond].append(0)
                continue
            
            if strand == "+":
                start_idx = region.index(exon_structure[2])
                extended_start_idx =  region.index(exon_structure[2] - add_us_ext[cond][cvg_idx])

                end_idx = region.index(exon_structure[3])
            else:
                start_idx = region.index(exon_structure[2])
                try:
                    end_idx = region.index(exon_structure[3])
                except ValueError:
                    end_idx = region.index(exon_structure[3] - 1) + 1
                try:
                    extended_end_idx = region.index(exon_structure[3] + add_us_ext[cond][cvg_idx])
                except ValueError:
                    # the end of the extended exon might not be accessible as position anymore
                    extended_end_idx = region.index(exon_structure[3] + add_us_ext[cond][cvg_idx] - 1) + 1

            if strand == "+":
                y = cvg[start_idx:end_idx]
                extended_y = cvg[extended_start_idx:end_idx]
            else:
                # for the mse calculation, the orientation of the coverage is irrelevant
                y = cvg[start_idx:end_idx]
                extended_y = cvg[start_idx:extended_end_idx]
            reg_mean = np.mean( y )
            reg_diff = y - reg_mean
            reg_mse = np.mean( reg_diff**2)
            mse_total_region[cond].append(reg_mse)
            
            # also calculate the mse 
            # with the exon upstream extension (might be used later)
            ext_reg_mean = np.mean( extended_y )
            ext_reg_diff = extended_y - ext_reg_mean
            ext_reg_mse = np.mean( ext_reg_diff**2)
            mse_total_ext_region[cond].append(ext_reg_mse)

    used_sites = find_sites_with_usage(polyA_sites,
                                       exon_structure[2],
                                       exon_structure[3],
                                       mse_total_region,
                                       mse_total_ext_region,
                                       bp_coverages_dict,
                                       valid_mean_cvg_dict,
                                       region,
                                       strand,
                                       options.min_coverage_region,
                                       options.best_break_point_upstream_extension, 
                                       options.mse_ratio_threshold, 
                                       add_us_ext)

    # append distal site
    used_sites.append(distal_site)

    ######
    # fivth part
    ######

    # go over the exon once more
    # check if any better global break point exists

    unsupported_break_point_cases = 0
    if len(used_sites) > 1:
        # create a list containing all genomic positions covered by the 
        # used poly(A) sites
        still_undetected = used_sites[:-1]

        # prepare individual dicts for every sample:

        # genomic pos' that are covered by sites in each sample
        used_sites_positions = {}
        # exon boundaries
        curr_ex_start = exon_structure[2]
        curr_ex_start_dict = {}
        for cond in bp_coverages_dict:
            curr_ex_start_dict[cond]=[]
            for i in range(len(bp_coverages_dict[cond])):
                # initialize the exon start with the annotated start pos
                # append the extended upstream region later
                curr_ex_start_dict[cond].append(curr_ex_start)
        curr_ex_end = exon_structure[3]
        curr_ex_end_dict = {}
        # same for exon end that is appended by the upstream region
        # later in the case of minus strand exons
        for cond in bp_coverages_dict:
            curr_ex_end_dict[cond]=[]
            for i in range(len(bp_coverages_dict[cond])):
                curr_ex_end_dict[cond].append(curr_ex_end)
        
        # treat the most proximal site separately;
        # here, for different samples the site might cover different genomic regions because of
        # variable upstream_extend_regions
        if strand == "+":
            start_pos = used_sites[0][1].start - options.best_break_point_upstream_extension
            end_pos = used_sites[0][1].end
            end_pos_idx = region.index(end_pos)
        else:
            start_pos = used_sites[0][1].start
            start_pos_idx = region.index(start_pos)
            end_pos = used_sites[0][1].end + options.best_break_point_upstream_extension
        for cond in bp_coverages_dict:
            # also, prepare the dict to store the covered genomic positions
            used_sites_positions[cond] = []
            
            for cvg_idx in range(len(bp_coverages_dict[cond])):
                used_sites_positions[cond].append({})

                if cond not in still_undetected[0][3]:
                    # this site is not supported by any sample of this condition
                    # hence, no dict preparation necessary
                    continue

                # only treat the current sample if it supports any of the
                # inferred poly(A) sites
                if not still_undetected[0][4][cond][cvg_idx]:
                    continue
                
                curr_add_us_ext = add_us_ext[cond][cvg_idx]
                if strand == "+":
                    if start_pos < curr_ex_start - curr_add_us_ext:
                        # start pos (site_start - best-break-point-extension is even more upstream
                        # than the extended exon; hence, use the extended exon as boundary
                        start_pos_idx = region.index( curr_ex_start - curr_add_us_ext )
                        curr_ex_start_dict[cond][cvg_idx] -= curr_add_us_ext
                    elif start_pos < curr_ex_start:
                        # start pos is more upstream than the unextended exon start;
                        # hence, use this start_pos but adjust exon start pos
                        start_pos_idx = region.index( start_pos )
                        curr_ex_start_dict[cond][cvg_idx] -= curr_add_us_ext
                    else:
                        # start_pos is within normal exon annotation
                        start_pos_idx = region.index(start_pos)
                else:
                    # same for negative strand exons; see comments above for explanation
                    if end_pos  >= curr_ex_end + curr_add_us_ext:
                        end_pos_idx = region.index( curr_ex_end + curr_add_us_ext -1 ) + 1
                        curr_ex_end_dict[cond][cvg_idx] += curr_add_us_ext
                    elif end_pos > curr_ex_end:
                        end_pos_idx = region.index( end_pos )
                        curr_ex_end_dict[cond][cvg_idx] += curr_add_us_ext
                    else:
                        end_pos_idx = region.index(end_pos)
                for i in range(start_pos_idx, end_pos_idx):
                    used_sites_positions[cond][cvg_idx][i] = 0
        
        # iterate over all other sites except the most proximal one
        # for these sites, the covered genomic positions are all the same for all samples
        for site_idx in range(1,len(still_undetected)):

            if strand == "+":
                start_pos = (used_sites[site_idx][1].start
                             - options.best_break_point_upstream_extension)

                # first: make the region as long as necessary but not longer than possible
                if start_pos < used_sites[site_idx - 1][1].end:
                    start_pos_idx = region.index( used_sites[site_idx - 1][1].end )
                else:
                    start_pos_idx = region.index(start_pos)

                end_pos = used_sites[site_idx][1].end
                end_pos_idx = region.index(end_pos)
            else:
                start_pos = used_sites[site_idx][1].start
                start_pos_idx = region.index(start_pos)
                
                end_pos = (used_sites[site_idx][1].end
                           + options.best_break_point_upstream_extension)

                # first: make the region as long as necessary but not longer than possible
                if end_pos > used_sites[site_idx - 1][1].start:
                    end_pos_idx = region.index( used_sites[site_idx - 1][1].start)
                else:
                    end_pos_idx = region.index(end_pos)

            # mark the genomic positions for all relevant samples
            for cond in bp_coverages_dict:
                if cond not in still_undetected[site_idx][3]:
                    # this site is not supported by any sample of this condition
                    # hence, no dict preparation necessary
                    continue
                for cvg_idx in range(len(bp_coverages_dict[cond])):
                    # only treat the current sample if it supports any of the
                    # inferred poly(A) sites
                    if not still_undetected[site_idx][4][cond][cvg_idx]:
                        continue
                    for i in range(start_pos_idx, end_pos_idx):
                        used_sites_positions[cond][cvg_idx][i] = site_idx

                        
        # if the proximal site is close to the exon start, the 
        # extension of the exon is used
        # during inference of the best break point
        # (this was ensured above when the exon
        # start positions were inferred individually for every sample)

        # iterate over all samples and infer the global best break point(s)
        # separately for each sample when this sample supports any
        # inferred poly(A) site
        for cond in bp_coverages_dict:
            for cvg_idx in range(len(bp_coverages_dict[cond])):
                sample_specific_still_undetected = []
                # collect the sites that are supported by the current sample
                for site_idx in range(len(still_undetected)):
                    if cond not in still_undetected[site_idx][3]:
                        # this site is not supported by any sample of this condition
                        sample_specific_still_undetected.append( None )
                        continue
                    if still_undetected[site_idx][4][cond][cvg_idx]:
                        sample_specific_still_undetected.append( still_undetected[site_idx] )
                    else:
                        sample_specific_still_undetected.append( None )

                if len([i for i in sample_specific_still_undetected if i is not None]) == 0:
                    # this sample does not support any poly(A) site;
                    # skip global break point inference for it
                    continue


                start_idx = region.index(curr_ex_start_dict[cond][cvg_idx])
                end_idx = region.index(curr_ex_end_dict[cond][cvg_idx] - 1) + 1

                global_mse_result = check_global_best_breaking_points(bp_coverages_dict[cond][cvg_idx],
                                                                      cond,
                                                                      cvg_idx,
                                                                      region,
                                                                      used_sites_positions[cond][cvg_idx],
                                                                      sample_specific_still_undetected,
                                                                      start_idx,
                                                                      end_idx,
                                                                      options.min_coverage_region)

                if isinstance(global_mse_result, tuple):
                    
                    if global_mse_result[1] is None:
                        syserr("[ERROR] No (further) break point found for exon: %s\n" 
                               % exon_id)
                        return ([], valid_mean_cvg_dict)
                    else:
                        # unsupported break point for current coverage detected
                        
                        # syserr("[ERROR] Unsupported break point %i in sample #%i (cond %s) for exon: %s. Exclude this sample for the current exon!\n" 
                        #        % (region[ global_mse_result[1] ], cvg_idx, cond, exon_id) )
                        # return ([], valid_mean_cvg_dict,add_us_ext)
                        valid_mean_cvg_dict[cond][cvg_idx] = -1
                        unsupported_break_point_cases += 1
                        # retract support for any site by this sample
                        for tmp_site_idx in range(len(used_sites) - 1):
                            if cond in used_sites[tmp_site_idx][3]:
                                used_sites[tmp_site_idx][4][cond][cvg_idx] = False
                        continue
    
    else:
        curr_ex_start = exon_structure[2]
        curr_ex_end = exon_structure[3]
        valid_break_point_dict = find_best_break_point_in_full_exon(bp_coverages_dict,
                                                                    valid_mean_cvg_dict,
                                                                    region,
                                                                    strand,
                                                                    curr_ex_start,
                                                                    curr_ex_end,
                                                                    options.min_coverage_region,
                                                                    options.mse_ratio_threshold,
                                                                    add_us_ext)

        # return distal site together with the valid_mean_cvg_dict!
        # NOT the valid_mean_cvg_dict!
        return (used_sites, valid_break_point_dict)

    if unsupported_break_point_cases > 0:
        syserr("[INFO] For exon %s, %i samples are not further processed due to unsupported global break points in their coverage\n" % (exon_id, unsupported_break_point_cases))
    return (used_sites, valid_mean_cvg_dict)

def get_distal_expression(used_sites,
                        exon_structure,
                        exon_id,
                        region,
                        bp_coverages_dict,
                        valid_break_point_found,
                        read_length,
                        min_coverage_region):
    '''Define the expression of the distal
    site if only this one is inferred
    '''

    expression_dict = {}

    strand = exon_structure[1]

    if exon_structure[3] - read_length - exon_structure[2] < min_coverage_region:
        syserr("[INFO] Single distal site of exon %s is very close to exon start. " + 
               "[No reasonable expression estimation possible. Skip exon\n" % str(exon_id))
        return None
    
    if strand == "+":
        start_idx = region.index( exon_structure[2] )
        end_idx = region.index( exon_structure[3] - read_length)
    else:
        end_idx = region.index( exon_structure[3] - 1) + 1
        start_idx = region.index( exon_structure[2]  + read_length)

    for cond in bp_coverages_dict:
        expression_dict[cond] = []
        for cvg_idx in range(len(bp_coverages_dict[cond])):

            # only infer an expression value if the sample
            # has valid coverage and no global break point
            # here, valid_break_point_found is a dict of True and False values
            # due to function find_best_break_point_in_full_exon
            # True: there is a valid break point!!!
            if valid_break_point_found[cond][cvg_idx]:
                expression_dict[cond].append(-1)
                continue
            
            cvg = bp_coverages_dict[cond][cvg_idx]
            mean_cvg = np.mean( cvg[start_idx:end_idx] )
            expr = float("{0:.2f}".format( mean_cvg ))
            expression_dict[cond].append(expr)

    return expression_dict
    
def get_relative_usages(used_sites,
                        exon_structure,
                        exon_id,
                        region,
                        bp_coverages_dict,
                        valid_mean_cvg_dict,
                        read_length,
                        min_coverage_region,
                        add_us_extension):
    """For the obtained sites, infer the relative 
    usage based on the mean expression of regions"""
    
    # define the indices that mark the segment boundaries
    # store them from proximal to distal
    strand = exon_structure[1]
    region_means = {}
    ds_corrected_means = {}
    relative_usages = {}
    expressions = {}

    for cond in bp_coverages_dict:
        region_means[cond] = []
        relative_usages[cond] = []
        expressions[cond] = []
        ds_corrected_means[cond] = []
        cvg_counter = -1
        for cvg in bp_coverages_dict[cond]:
            cvg_counter += 1

            if strand == "+":
                start_idx = region.index( exon_structure[2] )
                if used_sites[0][1].start - read_length - min_coverage_region < exon_structure[2]:
                    start_idx = region.index( exon_structure[2] - add_us_extension[cond][cvg_counter] )
                try:
                    end_idx = region.index( used_sites[0][1].start - read_length )
                except ValueError:
                    end_idx = region.index( exon_structure[2] - add_us_extension[cond][cvg_counter] + 1)

                if end_idx - start_idx < min_coverage_region:
                    syserr(("[INFO] Distance from proximal site to exon start of exon %s cond %s sample #%i " +
                            "is very small despite upstream extension (%i nt). For relative usage inference, " +
                            "the read length is not excluded from the upstream region of the proximal site\n")
                           % (exon_id, cond, cvg_counter, add_us_extension[cond][cvg_counter]))
                    end_idx = region.index( used_sites[0][1].start )
            else:
                if used_sites[0][1].end  + read_length + min_coverage_region > exon_structure[3]:
                    end_idx = region.index( exon_structure[3] + add_us_extension[cond][cvg_counter] - 1) + 1
                else:
                    end_idx = region.index( exon_structure[3] - 1) + 1

                try:
                    start_idx = region.index( used_sites[0][1].end  + read_length)
                except ValueError:
                    start_idx = end_idx = region.index( exon_structure[3] + add_us_extension[cond][cvg_counter] - 2) + 1

                if end_idx - start_idx < min_coverage_region:
                    syserr(("[INFO] Distance from proximal site to exon start of exon %s cond %s sample #%i " +
                            "is very small despite upstream extension (%i nt). For relative usage inference, " +
                            "the read length is not excluded from the upstream region of the proximal site\n")
                           % (exon_id, cond, cvg_counter, add_us_extension[cond][cvg_counter]))
                    start_idx = region.index( used_sites[0][1].end )
       
            region_means[cond].append([])
            relative_usages[cond].append([])
            expressions[cond].append([])
            ds_corrected_means[cond].append([])
            if (valid_mean_cvg_dict[cond][cvg_counter] == 0 or
                valid_mean_cvg_dict[cond][cvg_counter] == -1):
                region_means[cond][cvg_counter].append(0)
            else:
                y = cvg[start_idx:end_idx]
                region_means[cond][cvg_counter].append(np.mean(y))
            
    # iterate over all poly(A) sites (from prox to dist) and get the mean 
    for site_idx in range(1, len(used_sites)):
        if strand == "+":
            start_idx = region.index( used_sites[site_idx-1][1].end )
            end_idx = region.index( used_sites[site_idx][1].start - read_length )
        else:
            start_idx = region.index( used_sites[site_idx][1].end  + read_length )
            end_idx = region.index( used_sites[site_idx-1][1].start )
        
        # consistency check: is region big enough
        if end_idx - start_idx < min_coverage_region:
            syserr(("[ERROR] Upstream region of site %s is too small " +
                    "for reliable mean estimation\n") % used_sites[site_idx][0])
            return {}, {}
    
        for cond in bp_coverages_dict:
            counter = -1
            for cvg in bp_coverages_dict[cond]:
                counter += 1
                if (valid_mean_cvg_dict[cond][counter] == 0 or
                    valid_mean_cvg_dict[cond][counter] == -1):
                    region_means[cond][counter].append(0)
                else:
                    y = cvg[start_idx:end_idx]
                    region_means[cond][counter].append(np.mean(y))
                
    
    # now, all mean values are obtained
    # next: calculate downstream corrected mean values
    # from distal to proximal site
    for cond in region_means:
        for sample_arr_idx in range(len(region_means[cond])):
            if (valid_mean_cvg_dict[cond][sample_arr_idx] == 0 or
                valid_mean_cvg_dict[cond][sample_arr_idx] == -1):
                continue
            ds_corrected_means[cond][sample_arr_idx].append( region_means[cond][sample_arr_idx][-1] )
            
    for site_idx in reversed(range(len(used_sites) -1)):
        for cond in region_means:
            for sample_arr_idx in range(len(region_means[cond])):
                if (valid_mean_cvg_dict[cond][sample_arr_idx] == 0 or
                    valid_mean_cvg_dict[cond][sample_arr_idx] == -1):
                    continue
                downstream_corrected_mean = region_means[cond][sample_arr_idx][site_idx] - sum(ds_corrected_means[cond][sample_arr_idx])
                if downstream_corrected_mean < 0:
                    downstream_corrected_mean = 0
                ds_corrected_means[cond][sample_arr_idx].insert(0, downstream_corrected_mean)

    # iterate once more over the mean values and calculate actual relative usage values
    for site_idx in range(len(used_sites)):
        for cond in region_means:
            for sample_arr_idx in range(len(region_means[cond])):
                if (valid_mean_cvg_dict[cond][sample_arr_idx] == 0 or
                    valid_mean_cvg_dict[cond][sample_arr_idx] == -1):
                    relative_usages[cond][sample_arr_idx].append(-1.0)
                    expressions[cond][sample_arr_idx].append(-1.0)
                else:
                    rel_use = float("{0:.2f}".format( float(ds_corrected_means[cond][sample_arr_idx][site_idx]) / sum(ds_corrected_means[cond][sample_arr_idx]) * 100 ) )
                    expr = float("{0:.2f}".format( float(ds_corrected_means[cond][sample_arr_idx][site_idx]) ))
                    relative_usages[cond][sample_arr_idx].append(rel_use)
                    expressions[cond][sample_arr_idx].append(expr)
                
    return relative_usages, expressions

def process_exon_wrapper_function( input_tuple ):
    """Wrapper function to conduct the 
    actual processing of each exon and its sites"""

    exon, exon_structure, region, bp_coverages_dict, merged_exons, polyA_sites, exon_ext_dict, options = list(input_tuple)

    (used_sites, valid_mean_cvg_dict) = get_pA_set(exon, 
                                                   exon_structure, 
                                                   region, 
                                                   bp_coverages_dict, 
                                                   merged_exons, 
                                                   polyA_sites,
                                                   exon_ext_dict,
                                                   options)

    # ATTENTION:
    # valid_mean_cvg_dict might be the valid_break_point_dict if only the distal site
    # was quantified
    if len(used_sites) > 1:
        mean_cvg_indicator_set = set()
        for cond in valid_mean_cvg_dict:
            for i in valid_mean_cvg_dict[cond]:
                mean_cvg_indicator_set.add(i)
        if len(mean_cvg_indicator_set) == 1 and list(mean_cvg_indicator_set)[0] == 0:
            # no sample had enough coverage
            syserr("[INFO] For no sample global best break points and inferred sites matched for exon %s. Skipped exon!\n"
                   % str(exon))
            return (exon, [],[],[])
        elif len(mean_cvg_indicator_set) == 1 and list(mean_cvg_indicator_set)[0] == -1:
            # no sample had proper global best break points
            syserr("[INFO] For no sample global best break points and inferred sites matched for exon %s. Skipped exon!\n"
                   % str(exon))
            return (exon, [],[],[])
        elif (len(mean_cvg_indicator_set) == 2 and
              0 in mean_cvg_indicator_set and
              -1 in mean_cvg_indicator_set):
            # either the samples had too little coverage or their global best break point did not match
            syserr("[INFO] All sample were discarded either due to global best break point inconsistencies or low coverage for exon %s. Skipped exon!\n"
                   % str(exon))
            return (exon, [],[],[])

        # iterate over all sites except the distal one
        # it might occur that a site lost support by all its supporting coverages
        # due to global break point inconsistencies. -> Then, remove this site
        to_del_idx = []
        for site_idx in range(len(used_sites) - 1):
            supporting_samples = sum([i for b in used_sites[site_idx][3] for i in used_sites[site_idx][4][b] ] )
            if 0 == supporting_samples:
                to_del_idx.append(site_idx)

        for site_idx in sorted(to_del_idx)[::-1]:
            del used_sites[ site_idx ]
        
        if len(to_del_idx) > 0:
            syserr("[INFO] Number of sites was reduced for exon %s. They had no support by any sample anymore due to global break point inconsistencies\n"
                   % str(exon))

        if len(used_sites) == 1:
            # after filtering bogus sites
            # only the distal site remained -> delete this one as well
            used_sites = []

    if len(used_sites) > 1:
        relative_usage_dict, expression_dict = get_relative_usages(used_sites,
                                                                   exon_structure,
                                                                   exon,
                                                                   region,
                                                                   bp_coverages_dict,
                                                                   valid_mean_cvg_dict,
                                                                   options.read_length,
                                                                   options.min_coverage_region,
                                                                   exon_ext_dict)


        if len(relative_usage_dict) > 0:
            rel_use_array = []
            expression_array = []
            for site_idx in range(len(used_sites)):
                conditions_counter = {}
                rel_use_array.append([])
                expression_array.append([])
                for file_idx in range(len(options.coverages)):
                    condition = options.conditions[file_idx]
                    if condition not in conditions_counter:
                        conditions_counter[condition] = -1
                    conditions_counter[condition] += 1
                    sample = conditions_counter[condition]
                    rel_use_array[-1].append( relative_usage_dict[condition][ sample ][site_idx])
                    expression_array[-1].append( expression_dict[condition][ sample ][site_idx])

            return (exon, used_sites, rel_use_array, expression_array)

        else:
            # no relative usages obtained
            syserr("[ERROR] Failed to infer relative usages for %s\n" 
                   % str(exon) )
            return (exon, [],[], [])
    elif len(used_sites) == 1:
        # not multiple sites inferred
        # obtain the expression for the distal site in all samples
        # without global break point
        expression_dict = get_distal_expression(used_sites,
                                                exon_structure,
                                                exon,
                                                region,
                                                bp_coverages_dict,
                                                valid_mean_cvg_dict,
                                                options.read_length,
                                                options.min_coverage_region)

        if expression_dict is None:
            return (exon, [],[], [])
        else:
            # check the case that all samples have a valid global break point
            total_sample_cnt = 0
            for cond in expression_dict:
                for cvg_idx in expression_dict[cond]:
                    total_sample_cnt += 1
            if sum([i for c in expression_dict for i in expression_dict[c]]) == -total_sample_cnt:
                   # no sites inferred
                   syserr("[INFO] Distal site skipped due to conflicting global break points in all samples in exon %s. Skipped exon!\n"
                          % str(exon))
                   return (exon, [],[], [])
            
            expression_array = []
            conditions_counter = {}
            expression_array.append([])
            for file_idx in range(len(options.coverages)):
                condition = options.conditions[file_idx]
                if condition not in conditions_counter:
                    conditions_counter[condition] = -1
                conditions_counter[condition] += 1
                sample = conditions_counter[condition]
                expression_array[-1].append( expression_dict[condition][ sample ])
            return (exon, used_sites, None, expression_array)

    else:
        # no sites inferred
        syserr("[INFO] No site was inferred for %s. Skipped exon!\n"
               % str(exon))
        return (exon, [],[], [])

def main(options):
    """Main logic of the script"""

    #
    # IMPORTANT HACK START!
    #
    # This version of the script which calculates expression of polyA sites
    # requires a dirty hack in order to be incorporated into
    # a bigger, modularised snakemake pipeline.
    # Since we need Snakemake to evaluate ALL samples as wildcards in the DAG
    # we need to provide names of ALL samples as command-line arguments here.
    # However, we want to infer pas expression only for those samples which pass quality control.
    # Therefore we provide here additionally a TSV table with the set of those samples
    # which passed the QC (created in a previous snakemake module but
    # NOT AVAILABLE at the time of the call).
    # Therefore here, at the very beginning of the script
    # we need to artificially filter command-line args provided.
    # ~ Maciek
    #
    # parse the samples table:
    with open(options.design_file,"r") as design_f:
        design_lines = design_f.read().splitlines() # read the filtered samples table
    column_index = design_lines[0].split("\t").index("sample") 
    filtered_samples = [line.split("\t")[column_index] for line in design_lines[1:]] # get filtered samples names
    OK_indices_list = []
    # get CLArgs positional indices of those filtered samples
    for s in range(len(options.sample_names)):
        if options.sample_names[s] in filtered_samples:
            OK_indices_list.append(s)
    # filter all relevant CLArgs:
    options.sample_names = [options.sample_names[i] for i in OK_indices_list]
    options.coverages = [options.coverages[i] for i in OK_indices_list]
    options.conditions = [options.conditions[i] for i in OK_indices_list]
    options.ex_extension_files = [options.ex_extension_files[i] for i in OK_indices_list]
    #
    # IMPORTANT HACK END!
    #

    # initialize the pool for multimapping
    pool = multiprocessing.Pool( options.proc )

    expression_output_file_name = options.expressions_out
    distal_sites_output_file_name = options.distal_sites
    
    ###
    # load poly(A) sites
    ###

    # store the terminal exons
    terminal_exons_dict = {}

    # store the exon upstream extension size per sample
    term_ex_us_ex_dict = {}
    
    with open(options.polyA) as f:
        for line in f:
            F = line.strip().split("\t")
            chrom = F[0]
            strand = F[5]
            if F[8] not in terminal_exons_dict:
                # save stats for this exon and an empty list for the poly(A) sites
                # adjust the start coordinate of the exon to match
                # HTSeq and bed conventions
                exon_stats = F[8].split(":")
                exon_start = int(exon_stats[3]) - 1
                exon_end = int(exon_stats[4])
                terminal_exons_dict[ F[8] ] = [chrom, strand, exon_start, exon_end, [] ]
            # for each poly(A) site store: id, HTSeq.GenomicInterval of its location, protocol support
            terminal_exons_dict[ F[8] ][4].append( [ F[3], HTSeq.GenomicInterval( chrom, int(F[1]), int(F[2]), strand ), F[4] ] )
            term_ex_us_ex_dict[ F[8] ] = {}
            
    # go through list of all terminal exons and separate two lists:
    # 1. exons with clusters that assemble in close proximity
    # 2. exons with distinct clusters
    # clusters are expected to be sorted proximal -> distal
    for exon in terminal_exons_dict:
        strand = terminal_exons_dict[exon][1]
        list_of_close_clusters = []
        close_cluster_idxs = set()
        polyA_sites = terminal_exons_dict[ exon ][4]
        for cl_idx in range(1, len(polyA_sites)):
            prox_site = polyA_sites[cl_idx -1]
            dist_site = polyA_sites[cl_idx]
            if ( (strand == "+" and (dist_site[1].start - prox_site[1].end) < options.min_cluster_distance) or
            (strand == "-" and (prox_site[1].start - dist_site[1].end) < options.min_cluster_distance) ):
                # note these two sites as interferring
                if (cl_idx - 1) in close_cluster_idxs:
                    # simply update the set of close clusters
                    close_cluster_idxs.update([cl_idx])
                else:
                    # remove all clusters that are in the set so far
                    # store them separately as one unit
                    # update the set of closely spaced sites
                    if len(close_cluster_idxs) > 0:
                        if len(close_cluster_idxs) < 2:
                            syserr("[ERROR] set of closely spaced sites has less than 2 entries (%s) for exon %s\n" % (str(close_cluster_idxs), exon))
                            sys.exit(2) 
                        list_of_close_clusters.append(sorted(close_cluster_idxs))
                        close_cluster_idxs = set()
                    close_cluster_idxs.update([(cl_idx-1), cl_idx])
        # empty set if necessary
        if len(close_cluster_idxs) > 0:
            list_of_close_clusters.append(sorted(close_cluster_idxs))
                        
        if len(list_of_close_clusters) > 0:
            new_pA_list = []
            last_processed_pA = 0
            for close_clusters in list_of_close_clusters:
                for idx in range(last_processed_pA,close_clusters[0]):
                    new_pA_list.append(polyA_sites[idx])
                tmp_list = []
                for close_cl in close_clusters:
                    tmp_list.append(polyA_sites[close_cl])
                new_pA_list.append(tmp_list)
                last_processed_pA = close_clusters[-1] + 1
            for idx in range(last_processed_pA, len(polyA_sites)):
                new_pA_list.append(polyA_sites[idx])
            terminal_exons_dict[exon][4] = new_pA_list

    ###
    # load the coverages
    # and conditions
    ###

    coverage_pkl_dict = {}
    for pkl_idx in range(len(options.coverages)):
        extension_file = None
        pkl_file = options.coverages[pkl_idx]
        
        sample_basename = pkl_file.replace(".pkl", "")
        sample_basename = os.path.basename(sample_basename)
        for tmp_file in options.ex_extension_files:
            if sample_basename in tmp_file:
                extension_file = tmp_file
                
        if options.conditions[pkl_idx] not in coverage_pkl_dict:
            coverage_pkl_dict[options.conditions[pkl_idx]] = []
            # prepare the dict for the exon upstream extensions
            for ex in term_ex_us_ex_dict:
                term_ex_us_ex_dict[ex][options.conditions[pkl_idx]] = []

        if not pkl_file.endswith("pkl"):
            syserr("[ERROR] No proper pkl file as input: %s\n" % pkl_file)
            sys.exit(2)

        with open(pkl_file, 'rb') as input:
            cvg = cPickle.load(input)
            coverage_pkl_dict[options.conditions[pkl_idx] ].append(cvg)

        number_of_files = len( coverage_pkl_dict[options.conditions[pkl_idx] ] )
        if extension_file is not None and os.path.isfile(extension_file):
            with open(extension_file, "r") as  ext_file:
                for line in ext_file:
                    F = line.rstrip().split("\t")
                    # append the extension for all exons found in the file
                    term_ex_us_ex_dict[F[0]][ options.conditions[pkl_idx] ].append( int(F[1]) )
        # if an exon was not in the extension file
        # store 0 as extension
        for ex in term_ex_us_ex_dict:
            if len( term_ex_us_ex_dict[ex][options.conditions[pkl_idx] ] ) == number_of_files - 1:
                term_ex_us_ex_dict[ex][ options.conditions[pkl_idx] ].append(0)
            elif len( term_ex_us_ex_dict[ex][options.conditions[pkl_idx] ] ) != number_of_files:
                syserr(("[ERROR] Number of entries for exon extensions of ex %s does " +
                        "not match current number of processed pickle files for cond %s\n")
                       % (ex, options.conditions[pkl_idx]))
                sys.exit(2)

    # Finished reading input
    syserr("[INFO] %s Finished reading input\n" % time_now())
    
    ###
    # start iterating over all exons:
    ###

    # stats
    nr_skipped_exons = 0
    nr_no_negative_slope_fitted = 0
    nr_too_little_cov = 0
    
    data_entries = []

    for exon in terminal_exons_dict:
        
        exon_structure = terminal_exons_dict[ exon ][0:4]
        chrom = exon_structure[0]
        strand = exon_structure[1]
        start = exon_structure[2]
        end = exon_structure[3]

        polyA_sites = terminal_exons_dict[ exon ][4]

        max_ex_upstream_ext = max( [i for k in term_ex_us_ex_dict[exon] for i in term_ex_us_ex_dict[exon][k] ] )
        
        # extend the exon 3' boundary if the most distal site 
        # end more downstream than the annotated exon end
        if strand == "+":
            if isinstance(polyA_sites[-1][0], list) and polyA_sites[-1][-1][1].end > end:
                end = polyA_sites[-1][-1][1].end
            if not isinstance(polyA_sites[-1][0], list) and polyA_sites[-1][1].end > end:
                end = polyA_sites[-1][1].end
                
        else:
            if isinstance(polyA_sites[-1][0], list) and polyA_sites[-1][-1][1].start < start:
                start = polyA_sites[-1][-1][1].start
            if not isinstance(polyA_sites[-1][0], list) and polyA_sites[-1][1].start < start:
                start = polyA_sites[-1][1].start
            
        ### merge the coverages of this exon (used for distal site definition)
        if strand == "+":
            region = range( (start - max_ex_upstream_ext ), (end + options.ds_reg_for_no_coverage) )
        else:
            region = range( (start - options.ds_reg_for_no_coverage), (end + max_ex_upstream_ext) )
        bp_coverages_dict = {}
        merged_exons = np.zeros(end - start + options.ds_reg_for_no_coverage + max_ex_upstream_ext)
        for cond in coverage_pkl_dict:
            bp_coverages_dict[cond] = []
            for cvg in coverage_pkl_dict[cond]:
            
                # extend the 3' end to be able to check the downstream region as well
                # extend the 5' start to allow handling of proximal site close to the exon start
                if strand == "+":
                    curr_cov_arr = np.array( list( cvg[
                        HTSeq.GenomicInterval( chrom, 
                                               (start - max_ex_upstream_ext), 
                                               (end + options.ds_reg_for_no_coverage), 
                                               strand)
                    ]
                    )
                    )
                    
                else:
                    # do NOT reverse coverage array because
                    # it needs to match the region list
                    curr_cov_arr = np.array( list( cvg[
                        HTSeq.GenomicInterval( chrom, 
                                               (start - options.ds_reg_for_no_coverage), 
                                               (end + max_ex_upstream_ext), 
                                               strand)
                    ]
                    )
                    )
                # test
                # normalize the region so that the maximum is at 100
                #curr_cov_arr = curr_cov_arr / max(curr_cov_arr) * 100
                bp_coverages_dict[cond].append(curr_cov_arr)
                merged_exons = merged_exons + curr_cov_arr

        data_entries.append((exon, exon_structure, region, bp_coverages_dict, merged_exons, polyA_sites, term_ex_us_ex_dict[exon], options))
        
    result_tuples = pool.map( process_exon_wrapper_function, data_entries)

    with open(expression_output_file_name, "w") as ofile:
        with open( distal_sites_output_file_name, "w") as distal_out:

            # print header lines for the distal site and the expression file
            distal_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                "chrom",
                "start",
                "end",
                "pas",
                "score",
                "strand",
                "polyAsite_exon_idx",
                "nr_polyAsites_on_exon",
                "exon",
                "gene",
                "\t".join(options.sample_names) ) )

            ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                "chrom",
                "start",
                "end",
                "pas",
                "score",
                "strand",
                "polyAsite_exon_idx",
                "nr_polyAsites_on_exon",
                "exon",
                "gene",
                "\t".join(options.sample_names) ) )

            for res_tu in result_tuples:
                exon = res_tu[0]
                used_sites = res_tu[1]
                rel_use_array = res_tu[2]
                expression_array = res_tu[3]
                if len(used_sites) == 0 and len(rel_use_array) == 0:
                    continue
                elif len(used_sites) == 0 or (rel_use_array is not None and len(rel_use_array) == 0):
                    syserr(("[ERROR] length of list with inferred sites: %i\n" +
                            "length of list with relative usages: %i\n")
                           % (len(used_sites), len(rel_use_array)))
                    syserr("[ERROR] Only the case of both lists have zero length " +
                       "is expected\n")
                    sys.exit(2)
                elif len(used_sites) == 1:
                    # sanity check
                    if rel_use_array is not None:
                        syserr("[ERROR] Only distal site was inferred but quantification failed\n")
                        sys.exit(2)
                    expression_vals = "\t".join(str(x) for x in expression_array[0])
                    total_number_of_sites = len(used_sites)
                    transcript = exon.split(":")[0]
                    site = used_sites[0]
                    distal_out.write("%s\t%i\t%i\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n" % (site[1].chrom,
                                                                                       site[1].start,
                                                                                       site[1].end,
                                                                                       site[0],
                                                                                       site[2],
                                                                                       site[1].strand,
                                                                                       1,
                                                                                       total_number_of_sites,
                                                                                       exon,
                                                                                       transcript,
                                                                                       expression_vals))
                else: 
                    for idx in range(len(rel_use_array)):
                        rel_use_array[idx] = "\t".join(str(x) for x in rel_use_array[idx])
                        expression_array[idx] = "\t".join(str(x) for x in expression_array[idx])
                    total_number_of_sites = len(used_sites)
                    transcript = exon.split(":")[0]
                    cnt = 0
                    for site_idx in range(len(used_sites)):
                        site = used_sites[site_idx]
                        cnt += 1
                        sysout("%s\t%i\t%i\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n" % (site[1].chrom,
                                                                             site[1].start,
                                                                             site[1].end,
                                                                             site[0],
                                                                             site[2],
                                                                             site[1].strand,
                                                                             cnt,
                                                                             total_number_of_sites,
                                                                             exon,
                                                                             transcript,
                                                                             rel_use_array[site_idx]))
    
                        ofile.write("%s\t%i\t%i\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n" % (site[1].chrom,
                                                                                      site[1].start,
                                                                                      site[1].end,
                                                                                      site[0],
                                                                                      site[2],
                                                                                      site[1].strand,
                                                                                      cnt,
                                                                                      total_number_of_sites,
                                                                                      exon,
                                                                                      transcript,
                                                                                      expression_array[site_idx]))

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

        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
