"""
Abstract Class to Access Sequence, Quality from Sam Records
"""
from abc import ABCMeta, abstractmethod
import sam_constants

class SamSequence:
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_name(self):
        """
        :return str:  Return read name as seen in sam file.
        """

    @abstractmethod
    def is_intersect_slice(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Returns whether the read intersects the slice coordinates.
        :param int slice_start_wrt_ref_1based:  1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based:  1-based slice end position with respect to reference
        :return:  True if the sequence hits the slice (even with deletions, insertions, or non bases).  False otherwise.
        :rtype: bool
        """

    @abstractmethod
    def get_slice_intersect_coord(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Gets the 1-based position with respect to the reference of the intersection
        between the sam record sequence and the slice.
        If there is no intersection, returns (None, None)
        :param int slice_start_wrt_ref_1based: 1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based: 1-based slice end position with respect to reference
        :return tuple (int, int): (intersection start, intersection end)
        """

    @abstractmethod
    def get_seq_qual(self, do_pad_wrt_ref=False, do_pad_wrt_slice=False, do_mask_low_qual=False, q_cutoff=10,
                     slice_start_wrt_ref_1based=None, slice_end_wrt_ref_1based=None, do_insert_wrt_ref=False,
                     do_mask_stop_codon=False, stats=None, pad_space_btn_segments="N"):
        """
        Gets sequence, quality, AlignStats for the sam record.
        ASSUMES:  that sam cigar strings do not allow insertions at beginning or end of alignment (i.e.  local alignment as opposed to global alignment).

        :param bool do_pad_wrt_ref: Pad sequence with gaps with respect to the reference.
                        If slice_start_wrt_ref and slice_end_wrt_ref are valid, then pads the slice with respect to the reference.
        :param do_pad_wrt_slice: Ignored if do_pad_wrt_ref=True.  If True, then pads with gaps with respect to the slice.
        :param bool do_mask_low_qual: Mask bases with quality < q_cutoff with "N"
        :param int q_cutoff:  quality cutoff
        :param int slice_start_wrt_ref_1based:  If None, then whole sequence returned.   Otherwise, the slice 1-based start position with respect to the reference.
        :param int slice_end_wrt_ref_1based:  If None, then whole sequence returned.  Otherwise the slice end 1-based position with respect to the reference.
        :param bool do_insert_wrt_ref: Include insertions with respect to the reference.
                        If slice_start_wrt_ref and slice_end_wrt_ref are valid, then only includes inserts inside the slice.
        :param bool do_mask_stop_codon:  If True, then masks stop codons with "NNN".  Assumes that the reference starts that the beginning of a codon.
                Performs the stop codon masking after low quality base masking.
        :param AlignStats stats:  keeps track of stats.  Only counts inserts and quality if you allow inserts and mask quality.
        :param str pad_space_btn_segments:  The character to use in between aligned segments of a read.
        :return tuple (str, str, AlignStats):  (sequence, quality, AlignStats)
        """

    @abstractmethod
    def get_read_start_wrt_ref(self):
        """
        Gets the 1-based start position with respect to the reference of the unclipped portion of the read.
        :return:
        """

    @abstractmethod
    def get_read_end_wrt_ref(self):
        """
        Gets the 1-based end position with respect to the reference of the unclipped portion of the read.
        :return:
        """

    @staticmethod
    def do_pad(seq, seq_start_wrt_ref, seq_end_wrt_ref, pad_start_wrt_ref, pad_end_wrt_ref, pad_char=sam_constants.SEQ_PAD_CHAR):
        """
        Left and Right Pads the sequence.
        :param seq:
        :param seq_start_wrt_ref:
        :param seq_end_wrt_ref:
        :param pad_end_wrt_ref:
        :param str pad_char:  By default, pads with "-"
        :return:
        """
        if not seq_start_wrt_ref or not seq_end_wrt_ref:
            padded_seq = pad_char * (pad_end_wrt_ref - pad_start_wrt_ref + 1)
        else:
            left_pad_len = seq_start_wrt_ref  - pad_start_wrt_ref
            right_pad_len = pad_end_wrt_ref - seq_end_wrt_ref
            padded_seq = (pad_char * left_pad_len) + seq + (pad_char * right_pad_len)
        return padded_seq

    @classmethod
    def __subclasshook__(cls, C):
        return NotImplemented
