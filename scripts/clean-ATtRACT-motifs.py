"""
##############################################################################
#
#   Filter ATtRACT database and save motifs in a MotEvo-compatible format
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 23-11-2021
#   LICENSE: Apache_2.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import re
import os
import numpy as np
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
        "--dbpath", dest="dbpath", required=True, help="Path to the database textfile.",
    )
    parser.add_argument(
        "--pwmpath",
        dest="pwmpath",
        required=True,
        help="Path to the database pwm file.",
    )
    parser.add_argument(
        "--outdir", dest="outdir", required=True, help="Path for the output directory.",
    )
    return parser


##############################################################################


# custom function adapted from the initial work for:
# https://github.com/gruber-sciencelab/SMEAGOL
def read_pms_from_file(file, value_col="probs", lengths=False, transpose=False):
    """Function to read position matrices from a fasta-like file in Attract format.
    Args:
        pm_file (str): file containing PMs
        value_col (str): name for column containing PM values
        lengths (bool): lengths are provided in the file
        transpose (bool): transpose the matrix
    Returns:
        pandas dataframe containing PMs
    """
    # Read file
    pms = list(open(file, "r"))
    pms = [x.strip().split("\t") for x in pms]

    # Get matrix start and end positions
    starts = np.where([x[0].startswith(">") for x in pms])[0]
    assert starts[0] == 0
    ends = np.append(starts[1:], len(pms))

    # Get matrix IDs and values
    pm_ids = [l[0].strip(">") for l in pms if l[0].startswith(">")]
    if lengths:
        lens = np.array([l[1] for l in pms if l[0].startswith(">")]).astype("int")
        assert np.all(lens == ends - starts - 1)
    pms = [pms[start + 1 : end] for start, end in zip(starts, ends)]
    if transpose:
        pms = [np.transpose(np.array(x).astype("float")) for x in pms]
    else:
        pms = [np.array(x).astype("float") for x in pms]

    # Make dataframe
    return pd.DataFrame({"Matrix_id": pm_ids, value_col: pms})


# custom function adapted from:
# https://github.com/gruber-sciencelab/SMEAGOL
def trim_ppm(probs, information_content):
    """Function to trim non-informative columns from ends of a PPM.
    Args:
        probs (np.array): array containing PPM probability values
        frac_threshold (float): threshold (0-1) to filter out non-informative columns.
    Returns:
        result (np.array): array containing trimmed PPM.
    """
    pos_ic = position_wise_ic(probs, axis=1)
    to_trim = pos_ic <= information_content
    positions = list(range(probs.shape[0]))
    assert len(to_trim) == len(positions)
    # Trim from start
    while to_trim[0]:
        positions = positions[1:]
        to_trim = to_trim[1:]
    # Trim from end
    while to_trim[-1]:
        positions = positions[:-1]
        to_trim = to_trim[:-1]
    result = probs[positions, :]
    return result


# custom function adapted from:
# https://github.com/gruber-sciencelab/SMEAGOL
def position_wise_ic(probs, axis=1):
    """Function to calculate information content of each column in a PPM.
    Args:
        probs (np.array): array containing PPM probability values
    Returns:
        result (np.array): information content of each column in probs.
    """
    position_wise_entropy = np.apply_along_axis(entropy, axis=axis, arr=probs)
    result = 2 - position_wise_entropy
    return result


# custom function adapted from:
# https://github.com/gruber-sciencelab/SMEAGOL
def entropy(probs):
    """Function to calculate entropy of a PPM or column of a PPM.
    Args:
        probs (np.array): Array containing probability values
    Returns:
        result (float): Entropy value
    """
    result = -np.sum(probs * np.log2(probs))
    return result


def main():
    """Main body of the script."""

    # read in the motif annotation file
    db = pd.read_csv(options.dbpath, sep="\t")

    # keep non-mutated hsa RBPs
    db = db[db["Organism"] == "Homo_sapiens"]
    db = db[db["Mutated"] == "no"]

    # Unify name for all SELEX-experiments
    is_selex = np.array(
        [re.search("SELEX", _) is not None for _ in db["Experiment_description"]]
    )
    db.loc[is_selex, "Experiment_description"] = "SELEX"

    # Filter unnecessary columns, drop some duplication alreday at this point!
    # i.e. records linked to the same gene w/ the same pwm and experiment type
    # (probably these will be kmer variants from attract)
    cols = ["Gene_name", "Gene_id", "Matrix_id", "Experiment_description"]
    db = db[cols].drop_duplicates().reset_index(drop=True)

    # read in the motif PWM file
    pwms = read_pms_from_file(options.pwmpath, value_col="probs", lengths=True)

    # re-normalize probabilities
    pwms.probs = [x / np.expand_dims(np.sum(x, axis=1), 1) for x in pwms.probs]

    # merge annotations to pwms
    pwms = pwms.merge(db, on="Matrix_id").reset_index(drop=True)

    # Cluster by identical PWMs and drop all assigned to more than one gene
    duplicated_ids = set([])
    for i, row_i in pwms.iterrows():
        for j, row_j in pwms.iterrows():
            if i != j:  # compare two distinct rows
                # if the matrices are equal
                if np.array_equal(pwms.at[i, "probs"], pwms.at[j, "probs"]):
                    # if they are assigned to a different gene
                    # - this is crap and we get rid of this PWM
                    if row_i["Gene_name"] != row_j["Gene_name"]:
                        duplicated_ids.add(row_i["Matrix_id"])
                        duplicated_ids.add(row_j["Matrix_id"])
                    # if they are assigned to the same gene
                    # - this is a duplication in the database and
                    # we want to keep just one record
                    else:
                        # mark the second occurance as duplication event
                        # and keep the first pwm in
                        duplicated_ids.add(row_j["Matrix_id"])
    unique = pwms[~pwms["Matrix_id"].isin(duplicated_ids)].copy()

    # Trim PWMs by position-wise IC
    unique["orig_len"] = [p.shape[0] for p in unique.probs]
    unique["trimmed_probs"] = unique.probs.apply(
        trim_ppm, args=(0.0,)
    )  # trim non-informative columns from ends of a PWM
    unique["trimmed_len"] = [p.shape[0] for p in unique.trimmed_probs]

    # Filter by length of the core of the motif (trimmed len)
    # This step essentially removes long motifs which do not have at leas 4len "core"
    # also discards "essentially kmers" - streches of 12+ polarized matrices
    lentrimmed_kept = unique[
        ((unique.trimmed_len >= 4) & (unique.trimmed_len <= 7))
    ].copy()
    lentrimmed_kept.drop(columns=["probs", "orig_len"], inplace=True)
    lentrimmed_kept.rename(
        columns={"trimmed_probs": "probs", "trimmed_len": "len"}, inplace=True
    )

    # Filter by entropy
    lentrimmed_kept["entropy"] = lentrimmed_kept.probs.apply(entropy)
    max_entropy = 10
    entropy_kept = lentrimmed_kept[(lentrimmed_kept.entropy <= max_entropy)].copy()

    # Some gene names contain spaces in the string!
    # this messes up the ID downstream
    # strip whitespaces!
    pwms = entropy_kept.reset_index(drop=True)
    pwms["Gene_name"] = pwms.apply(lambda row: row["Gene_name"].strip(), axis=1)

    # save motifs in separate textfiles
    os.makedirs(options.outdir)
    for i, row in pwms.iterrows():
        fullname = row["Gene_name"] + "_" + row["Matrix_id"]
        with open(os.path.join(options.outdir, fullname), "w") as f:
            f.write("//" + os.linesep)
            f.write("NA " + fullname + os.linesep)
            f.write("\tA\tC\tG\tT" + os.linesep)
            P = row["probs"] * 100
            for i in range(1, row["len"] + 1):
                stri = str(i)
                if len(stri) == 1:
                    stri = "0" + stri
                f.write(stri)
                f.write(
                    "\t"
                    + str(round(P[i - 1, 0], 3))
                    + "\t"
                    + str(round(P[i - 1, 1], 3))
                    + "\t"
                    + str(round(P[i - 1, 2], 3))
                    + "\t"
                    + str(round(P[i - 1, 3], 3))
                    + os.linesep
                )
            f.write("//" + os.linesep)


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
