# Finds duplicate reads
# Before slicing, find all duplicates across entire reference
# Reads in query sort order (easier - no change to pipeline)
# Merge reads, quality mask, conflict mask.
# Keep track of coordinates of merged reads and sequence  --> This will be massive in memory for large sized organisms.
# Use a dict for now because easiest.
# Key = tuple (left-most coordinate that hits reference, sequence with external gaps trimmed)
# Don't both keeping rightmost coordinate, since calling len(str) is O(1) in python :  http://stackoverflow.com/questions/1115313/cost-of-len-function
# In order to be duplicate, reads must have exact same sequence (with external gaps trimmed) and same start coordinate
# For every window slice, instead of iterating through sam, we iterate through our dict of merged reads.
# Check if start coordinate of merged read is on or after slice start.
# Calculate the end coordinate of merged read is before or on slice end
# If merged read is in slice then put it into window
# This means that we will have a separate copy of our merged read dict for every parallel window.  Sometimes we have sam files that are 300MB big just for 1 patient sample for HIV.  This will probably expand to 300MB x 10 = 3GB for each window.  This will limit the number of parallel windows we can execute on 1 node.
#

import sam_handler
from collections import namedtuple
import logging

LOGGER = logging.getLogger(__name__)

# start = 1 based start position with respect to reference
# seq = merged sequence, external gaps trimmed
UniqSeq = namedtuple("UniqSeq", ["start", "seq"], verbose=(LOGGER.level == logging.DEBUG))
# read = readname
# score = sum of quality scores for bases that align to reference and aren't masked for low quality or conflict
# sam_idx = 0-based index of mate1 record in sam
UniqSeqScore = namedtuple("UniqSeqScore", ["read", "score", "sam_idx"], verbose=(LOGGER.level == logging.DEBUG))


def find_dups(sam_filename, ref, out_csv_filename, mapping_cutoff, read_qual_cutoff, max_prop_N, breadth_thresh, is_insert=False):
    """
    Goes through the sam, merges paired records, and finds merged sequences that are exact duplicates of other merged sequences.
    To be exact duplicate, both merged sequences must have the same start coordinates and matching sequence
    (including same gaps, N's, and bases).

    Keeps track of all unique sequences encountered in a dict.
    Key = tuple (left coordinate wrt ref, sequence with external gaps trimmed)
    This will be a memory hog.

    :param sam_filename:
    :param ref:
    :param str out_csv_filename:  output file
    :param mapping_cutoff:
    :param read_qual_cutoff:
    :param max_prop_N:
    :param breadth_thresh:
    :param is_insert:
    :return {UniqSeq: :  dict
    """


    sam_handler.create_msa_slice_from_sam(sam_filename=sam_filename, ref=ref, , mapping_cutoff, read_qual_cutoff, max_prop_N,
                                  breadth_thresh, start_pos=0, end_pos=0, is_insert=False, is_mask_stop_codon=False,
                                  ref_len=0)