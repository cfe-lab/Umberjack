import logging
import sys
import sam_constants
import sam.align_stats
import math
from sam.sam_record import SamRecord as SamRecord

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

class PairedRecord:
    """
    Only for mates that hit the same reference.
    """

    def __init__(self, mate1_record, mate2_record):
        self.mate1 = mate1_record
        self.mate2 = mate2_record
        if self.mate1.rname != self.mate2.rname:
            raise ValueError("mate1_record and mate2_record  must align to the same reference")
        self.mate1.fill_mate(self.mate2)
        self.read_start_wrt_ref = None
        self.read_end_wrt_ref = None


    def get_name(self):
        return self.mate1.qname

    def get_read_start_wrt_ref(self):
        """
        Gets the 1-based position start of the merged read.
        :return:
        """
        if not self.read_start_wrt_ref:
            self.read_start_wrt_ref = min(self.mate1.get_seq_start_wrt_ref(), self.mate2.get_seq_start_wrt_ref())
        return self.read_start_wrt_ref


    def get_read_end_wrt_ref(self):
        """
        Gets the 1-based position end of the merged read.
        :return:
        """
        if not self.read_end_wrt_ref:
            self.read_end_wrt_ref = max(self.mate1.get_seq_end_wrt_ref(), self.mate2.get_seq_end_wrt_ref())
        return self.read_end_wrt_ref


    def is_slice_intersect_read(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Returns whether the slice intersects the read
        :param slice_start_wrt_ref_1based:
        :param slice_end_wrt_ref_1based:
        :return:
        """
        return slice_start_wrt_ref_1based <= self.get_read_end_wrt_ref() and slice_end_wrt_ref_1based >= self.get_read_start_wrt_ref()


    def is_between_mates(self, pos_wrt_ref_1based):
        """
        Returns whether mates have a gap in between them and the position is inside that gap.
        :param pos_wrt_ref_1based:  1-based position with respect to reference
        :return:
        """
        return ((self.mate1.get_seq_end_wrt_ref() < pos_wrt_ref_1based < self.mate2.get_seq_start_wrt_ref()) or
                (self.mate2.get_seq_end_wrt_ref() < pos_wrt_ref_1based < self.mate1.get_seq_start_wrt_ref()))


    def get_read_slice_intersect(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Gets the 1-based position with respect to the reference of the intersection between the read and slice.
        If there is no intersection, returns (None, None)
        :param slice_start_wrt_ref_1based:
        :param slice_start_wrt_ref_1based:
        :return tuple (int, int): (intersection start, intersection end)
        """
        read_slice_intersect_start_wrt_ref = None
        read_slice_intersect_end_wrt_ref = None
        if self.is_slice_intersect_read(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
            # 1-based position with respect to reference of the start of intersection of the read fragment and slice
            read_slice_intersect_start_wrt_ref = max(self.get_read_start_wrt_ref(), slice_start_wrt_ref_1based)
            # 1-based position with respect to reference of the end of intersection of read fragment and slice
            read_slice_intersect_end_wrt_ref = min(self.get_read_end_wrt_ref(), slice_end_wrt_ref_1based)
        return read_slice_intersect_start_wrt_ref, read_slice_intersect_end_wrt_ref



    def is_in_mate_overlap(self, pos_wrt_ref_1based):
        """
        Returns whether the 1-based position with respect to the reference is region of mate overlap.
        :param pos_wrt_ref_1based:  1-based position with respect to the reference
        :return:
        """
        return ((self.mate1.get_seq_start_wrt_ref() <= pos_wrt_ref_1based <= self.mate1.get_seq_end_wrt_ref()) and
                (self.mate2.get_seq_start_wrt_ref() <= pos_wrt_ref_1based <= self.mate2.get_seq_end_wrt_ref()))

    @staticmethod
    def is_in_slice(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based, pos_wrt_ref_1based):
        """
        Returns whether the position is within the slice
        :param slice_start_wrt_ref_1based:
        :param slice_end_wrt_ref_1based:
        :param pos_wrt_ref_1based:
        :return:
        """
        return slice_start_wrt_ref_1based <= pos_wrt_ref_1based <= slice_end_wrt_ref_1based


    @staticmethod
    def calc_q_cutoff_overlap(q_cutoff):
        # If overlapping mate bases agree, either mate quality can be lower as long as the
        #   Probability of single base error at the quality cutoff >= Probability of both mates having same erroneous base
        # Quality = -10 log P(error),  P(error) = 1e-(Quality/10)  by Phred Score definition.
        # Assume that the probability of a wrong base is independent across mates, positions.
        # Say the correct base = G, the permutations of wrong overlapping bases = {AA, AC, AT, CA, CC, CT, TA, TC, TT}.
        # ==> There are 3 permutations where the overlapping bases are the same; There are 9 possible permutations of wrong overlapping bases.
        #
        # P(single base error at q_cutoff) >= P(both mates same base & both mates wrong base)
        # P(single base error at q_cutoff) >= P(both mates wrong base) P(both mates same base | both mates wrong base)
        # 1e-(q_cutoff/10) >= P(mate1 1-base error) P(mate2 1-base error) [3 matching overlaps / 9 possible permutations of overlapping wrong bases)
        # 1e-(q_cutoff/10) >= 1e-q1/10 *  1e-q2/10  * 1/3
        # 3e-(q_cutoff/10) >= 1e-q1/10 *  1e-q2/10
        # log3 - q_cutoff/10 >= -q1/10 - q2/10
        # -10log3 + qcutoff <= q1 + q2
        # qcutoff_overlap = -10log3 + qcutoff <= q1 + q2
        q_cutoff_overlap = -10 * math.log10(3) + q_cutoff
        return q_cutoff_overlap


    def get_read_inserts(self, slice_start_wrt_ref_1based=None, slice_end_wrt_ref_1based=None):
        """
        Combines the dicts for mate1 and mate2 inserts into a single dict.  Ignores inserts outside of the slice.
        :param slice_start_wrt_ref_1based:
        :param slice_end_wrt_ref_1based:
        :return:
        """
        merge_rpos_to_insert = dict()  # {pos : { 1: (seq1, qual1), 2: (seq1, qual2)}
        for mate, insert_dict in [(1, self.mate1.get_insert_dict()), (2, self.mate2.get_insert_dict())]:
            for insert_pos, (insert_seq, insert_qual) in insert_dict.iteritems():
                # Ignore inserts outside of slice
                if PairedRecord.is_in_slice(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based-1, insert_pos):
                    if not merge_rpos_to_insert.get(insert_pos):
                        merge_rpos_to_insert[insert_pos] = dict()
                    if not merge_rpos_to_insert[insert_pos].get(mate):
                        merge_rpos_to_insert[insert_pos][mate] = dict()
                    merge_rpos_to_insert[insert_pos][mate] = (insert_seq, insert_qual)

        return merge_rpos_to_insert


    def merge_inserts(self, sliced_mseq, q_cutoff, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based,
                           stats=None):
        """
        Assumes that sliced_mseq only contains the sliced portion of the merged read sequence.

        :param sliced_mseq:
        :param q_cutoff:
        :param slice_start_wrt_ref_1based:
        :param slice_end_wrt_ref_1based:
        :return:
        """
        # track stats about inserts
        if not stats:
            stats = sam.align_stats.AlignStats()

        q_cutoff_overlap = PairedRecord.calc_q_cutoff_overlap(q_cutoff)

        read_insert_dict = self.get_read_inserts(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based)
        merge_rpos_to_insert = dict()

        # For the positions in which the mates overlap, check if the inserts are the same in both mates
        # insert_pos1_wrt_ref_1based is the 1-based position wrt ref before after the insertion
        for insert_pos_wrt_ref_1based, mate_dict in read_insert_dict.iteritems():
            stats.total_insert_blocks += 1
            # Does this insert fall into the mate overlap?
            if self.is_in_mate_overlap(insert_pos_wrt_ref_1based):
                insert_seq1, insert_qual1 = mate_dict.get(1, ("", ""))
                insert_seq2, insert_qual2 = mate_dict.get(2, ("", ""))

                # If the inserted position exists in one mate, then remove it.
                # If the inserted position exists in both mates, but they disagree on the base, take the higher quality base > q_cutoff.
                # If the inserted position exists in both mates, but they disagree on the base, mask it if both low quality.
                masked_insert_seq = ""
                length = max(len(insert_seq1), len(insert_seq2))
                for i in range(0, length):
                    stats.total_inserts += 1
                    ichar1 = None
                    ichar2 = None
                    iqual1 = None
                    iqual2 = None

                    if i < len(insert_seq1):
                        ichar1 = insert_seq1[i]
                        iqual1 = ord(insert_qual1[i])-sam_constants.PHRED_SANGER_OFFSET

                    if i < len(insert_seq2):
                        ichar2 = insert_seq2[i]
                        iqual2 = ord(insert_qual2[i])-sam_constants.PHRED_SANGER_OFFSET

                    if ichar1 == ichar2:
                        stats.total_insert_agree += 1
                        if iqual1 + iqual2 >= q_cutoff_overlap:
                            masked_insert_seq += ichar1
                            stats.total_insert_agree_hi_qual += 1
                        else:
                            masked_insert_seq += "N"
                            stats.total_insert_agree_lo_qual += 1
                    elif ichar2 and not ichar1 :  # Only mate2 has insert here.  Exclude insert.
                        stats.total_insert_conflict += 1
                        if iqual2 >= q_cutoff:
                            stats.total_insert_conflict_hino_qual += 1
                        else:
                            stats.total_insert_conflict_lono_qual += 1
                    elif ichar1 and not ichar2:  # Only mate1 has insert here.  Exclude insert.
                        stats.total_insert_conflict += 1
                        if iqual1 >= q_cutoff:
                            stats.total_insert_conflict_hino_qual += 1
                        else:
                            stats.total_insert_conflict_lono_qual += 1
                    elif ichar2 != ichar1:  # Mate1 and Mate2 both have inserts here, but bases conflict.
                        stats.total_insert_conflict += 1

                        if iqual1 < q_cutoff and iqual2 < q_cutoff:  # Both inserts low quality, mask it
                            masked_insert_seq += "N"
                            stats.total_insert_conflict_lo_qual += 1
                        elif iqual1 > iqual2 >= q_cutoff:  # both inserts high quality, but insert1 higher
                            masked_insert_seq += ichar1
                            stats.total_insert_conflict_hi_qual += 1
                        elif iqual2 > iqual1 >= q_cutoff:  # both sequences high quality, but seq2 higher
                            masked_insert_seq += ichar2
                            stats.total_insert_conflict_hi_qual += 1
                        elif iqual1 >= q_cutoff > iqual2:  # seq1 high quality, seq2 low quality
                            masked_insert_seq += ichar1
                            stats.total_insert_conflict_hilo_qual += 1
                        elif iqual2 >= q_cutoff > iqual1:  # seq2 high quality, seq2 low quality
                            masked_insert_seq += ichar2
                            stats.total_insert_conflict_hilo_qual += 1
                        elif iqual1 == iqual2 >= q_cutoff:  # both seq high equal quality, mask it
                            masked_insert_seq += "N"
                            stats.total_insert_conflict_equal_hi_qual += 1

                merge_rpos_to_insert[insert_pos_wrt_ref_1based] = masked_insert_seq
            else:
                # If the the insert is in a region that does not overlap with the mate, exclude it if it's low quality
                masked_insert_seq = ""
                if mate_dict.get(1):
                    insert_seq, insert_qual = mate_dict.get(1)
                else:
                    insert_seq, insert_qual = mate_dict.get(2)
                for i, ichar in enumerate(insert_seq):
                    stats.total_insert_1mate += 1
                    stats.total_inserts += 1
                    iqual = ord(insert_qual[i])-sam_constants.PHRED_SANGER_OFFSET
                    if iqual >= q_cutoff:
                        masked_insert_seq += ichar
                        stats.total_insert_1mate_hi_qual += 1
                    else:
                        stats.total_insert_1mate_lo_qual += 1

                merge_rpos_to_insert[insert_pos_wrt_ref_1based] = masked_insert_seq




        mseq_with_inserts = ""
        last_insert_pos_0based_wrt_mseq = -1  # 0-based position wrt result_seq before the previous insertion
        # insert_pos_wrt_ref: 1-based reference position before the insertion
        for insert_1based_pos_wrt_ref, insert_seq in merge_rpos_to_insert.iteritems():
            # 0-based position wrt sliced_mseq right before the insertion
            insert_pos_0based_wrt_mseq =  insert_1based_pos_wrt_ref - slice_start_wrt_ref_1based
            mseq_with_inserts += sliced_mseq[last_insert_pos_0based_wrt_mseq+1:insert_pos_0based_wrt_mseq+1] + insert_seq
            last_insert_pos_0based_wrt_mseq = insert_pos_0based_wrt_mseq

        # TODO:  assume that bowtie will not let inserts at beginning/end of align, but check
        mseq_with_inserts += sliced_mseq[last_insert_pos_0based_wrt_mseq+1:len(sliced_mseq)]

        if merge_rpos_to_insert:
            LOGGER.debug("qname=" + self.mate1.qname + " insert stats:\n" + stats.dump_insert_stats())

        return mseq_with_inserts, stats



    def merge_match(self, q_cutoff=10, pad_space_btn_mate="N", slice_start_wrt_ref_1based=None, slice_end_wrt_ref_1based=None, stats=None):
        """
        Merge bases in reads that align as a match to the reference (as opposed to bases that align as an insert to the ref).
        :param q_cutoff:
        :param slice_end_wrt_ref_1based:  If None or slice_start_wrt_ref_1based is None, then uses the read end wrt ref
        :param slice_start_wrt_ref_1based: If None or slice_end_wrt_ref_1based is None, then uses the read start wrt ref
        :param stats:
        :return:
        """
        mseq = ""

        if not slice_start_wrt_ref_1based or not slice_end_wrt_ref_1based:
            slice_start_wrt_ref_1based = self.read_start_wrt_ref()
            slice_end_wrt_ref_1based = self.read_end_wrt_ref()

        if slice_start_wrt_ref_1based > slice_end_wrt_ref_1based:
            raise ValueError("Slice start should be <= slice end")

        # Extract portion of reads that fit within the intersection of the read fragment & slice.
        # Pad the extracted sequences just enough so that they line up within the intersection.
        # Do not pad with respect to the reference yet.  Do not pad with respect to the slice yet.
        # Pop out insertions with respect to the reference so that it is easier to find the sequence coordinates wrt ref.
        # Keep track of the insertions and the ref pos right before the insertion.
        seq1, qual1, stats = self.mate1.get_seq_qual(do_pad_wrt_ref=False,
                                                       do_pad_wrt_slice=True,
                                                       do_mask_low_qual=False,
                                                       q_cutoff=q_cutoff,
                                                       do_insert_wrt_ref=False,
                                                       slice_start_wrt_ref_1based=slice_start_wrt_ref_1based,
                                                       slice_end_wrt_ref_1based=slice_end_wrt_ref_1based,
                                                       stats=stats)
        seq2, qual2, stats = self.mate2.get_seq_qual(do_pad_wrt_ref=False,
                                                       do_pad_wrt_slice=True,
                                                       do_mask_low_qual=False,
                                                       q_cutoff=q_cutoff,
                                                       do_insert_wrt_ref=False,
                                                       slice_start_wrt_ref_1based=slice_start_wrt_ref_1based,
                                                       slice_end_wrt_ref_1based=slice_end_wrt_ref_1based,
                                                       stats=stats)

        if len(seq1) != len(seq2) or len(qual1) != len(qual2) or len(qual1) != len(seq1):
            raise ValueError("Expect left and right pad mate sequences and qualities such that they line up. " +
                                " mate1=" + self.mate1.qname + " mate2=" + self.mate2.qname +
                                " len_seq1={} len_seq2={},  len_qual1={} len_qual2={}".format(len(seq1),
                                                                                              len(seq2),
                                                                                              len(qual1),
                                                                                              len(qual2)))

        # modified cutoff for overlapping matching bases.  Due to the match, this is lower than the q_cutoff.
        q_cutoff_overlap = PairedRecord.calc_q_cutoff_overlap(q_cutoff)

        for i in range(0, len(seq2)):
            pos_wrt_ref = i + slice_start_wrt_ref_1based
            q1 = ord(qual1[i])-sam_constants.PHRED_SANGER_OFFSET
            q2 = ord(qual2[i])-sam_constants.PHRED_SANGER_OFFSET


            # Both mates have gap at this position.
            if seq1[i] == seq2[i] == sam_constants.SEQ_PAD_CHAR:
                # Change to special pad character if the position is inside the read fragment
                # so that multiple sequence alignment realizes that there is a base here as instead of a gap wrt ref.
                if self.is_between_mates(pos_wrt_ref):
                    mseq += pad_space_btn_mate
                else:
                    mseq += sam_constants.SEQ_PAD_CHAR

            # only mate2 has real base here
            elif seq1[i] == sam_constants.SEQ_PAD_CHAR and seq2[i] != sam_constants.SEQ_PAD_CHAR:
                stats.total_match_1mate += 1
                if q2 >= q_cutoff:
                    mseq += seq2[i]
                    stats.total_match_1mate_hi_qual += 1
                else:
                    mseq += "N"
                    stats.total_match_1mate_lo_qual += 1

            # only mate1 has real base here
            elif seq2[i] == sam_constants.SEQ_PAD_CHAR and seq1[i] != sam_constants.SEQ_PAD_CHAR:
                stats.total_match_1mate += 1
                if q1 >= q_cutoff:
                    mseq += seq1[i]
                    stats.total_match_1mate_hi_qual += 1
                else:
                    mseq += "N"
                    stats.total_match_1mate_lo_qual += 1

            # Both mates have the same base at this position
            elif seq1[i] == seq2[i] and seq1[i] != sam_constants.SEQ_PAD_CHAR:
                stats.total_match_nonconflict += 1
                if q_cutoff_overlap <= q1 + q2:
                    mseq += seq1[i]
                    stats.total_match_nonconflict_hi_qual += 1
                else:
                    mseq += "N"
                    stats.total_match_nonconflict_lo_qual += 1

            # Both mates have disagreeing bases at this position: take the high confidence
            elif seq1[i] != seq2[i]:
                stats.total_match_conflict += 1
                if q1 < q_cutoff and q2 < q_cutoff:  # both sequences low quality
                    mseq += "N"
                    stats.total_match_conflict_lo_qual += 1
                elif q1 > q2 >= q_cutoff:  # both sequences high quality, but seq1 higher
                    mseq += seq1[i]
                    stats.total_match_conflict_hi_qual += 1
                elif q1 >= q_cutoff > q2:  # seq1 high quality, seq2 low quality
                    mseq += seq1[i]
                    stats.total_match_conflict_hilo_qual += 1
                elif q2 > q1 >= q_cutoff:  # both sequences high quality, but seq2 higher
                    mseq += seq2[i]
                    stats.total_match_conflict_hi_qual += 1
                elif q2 >= q_cutoff > q1:  # seq2 high quality, seq2 low quality
                    mseq += seq2[i]
                    stats.total_match_conflict_hilo_qual += 1
                elif q1 == q2 >= q_cutoff:
                    mseq += "N"
                    stats.total_match_conflict_equal_hi_qual += 1
                else:
                    raise ValueError("We should never get here.  Unanticipated use case for merging sam reads")

            else:
                raise ValueError("We should never get here.  Unanticipated use case for merging sam reads")

        return mseq


    def merge_sam_reads(self, q_cutoff=10, pad_space_btn_mates="N", do_insert_wrt_ref=False, do_pad_wrt_ref=True,
                        do_pad_wrt_slice=False, slice_end_wrt_ref_1based=None, slice_start_wrt_ref_1based=None,
                        stats=None):
        """
        Merge two sequences that overlap over some portion (paired-end
        reads).  Using the positional information in the SAM file, we will
        know where the sequences lie relative to one another.  In the case
        that the base in one read has no complement in the other read
        (in partial overlap region), take that base at face value.

        Allows insertions.

        :param pad_space_btn_mates:
        :param do_pad_wrt_slice:
        :return:  merged paired-end read
        :rtype : str
        :param int q_cutoff: quality cutoff below which a base is converted to N if there is no consensus between the mates.
        :param bool do_insert_wrt_ref: whether insertions with respect to reference is allowed.
        """

        mseq = ''
        if not stats:
            stats = sam.align_stats.AlignStats()

        if not slice_start_wrt_ref_1based:
            slice_start_wrt_ref_1based = 1
        if not slice_end_wrt_ref_1based:
            slice_end_wrt_ref_1based = self.mate1.ref_len

        # If the slice does not intersect either mate,
        # Then just return empty string or padded gaps wrt ref or slice as desired.
        if not self.is_slice_intersect_read(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
            if do_pad_wrt_ref:
                mseq = SamRecord.do_pad(mseq, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                        pad_start_wrt_ref=1, pad_end_wrt_ref=self.mate1.ref_len,
                                        pad_char=sam_constants.SEQ_PAD_CHAR)
            elif do_pad_wrt_slice:
                mseq = SamRecord.do_pad(mseq, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                        pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                        pad_char=sam_constants.SEQ_PAD_CHAR)
            return mseq, stats


        # 1-based position with respect to reference of intersection between the read fragment and slice
        read_slice_xsect_start_wrt_ref, read_slice_xsect_end_wrt_ref = self.get_read_slice_intersect(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based)

        mseq = self.merge_match(q_cutoff=10, pad_space_btn_mate=pad_space_btn_mates,
                         slice_start_wrt_ref_1based=read_slice_xsect_start_wrt_ref,
                         slice_end_wrt_ref_1based=read_slice_xsect_end_wrt_ref, stats=stats)


        if do_insert_wrt_ref:
            mseq, stats = self.merge_inserts(sliced_mseq=mseq, q_cutoff=q_cutoff,
                                                  slice_start_wrt_ref_1based=read_slice_xsect_start_wrt_ref,
                                                  slice_end_wrt_ref_1based=read_slice_xsect_end_wrt_ref, stats=stats)
        # Now pad with respect to reference
        if do_pad_wrt_ref:
            mseq = SamRecord.do_pad(seq=mseq, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                    seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                    pad_start_wrt_ref=1, pad_end_wrt_ref=self.mate1.ref_len)

        elif do_pad_wrt_slice:
            mseq = SamRecord.do_pad(seq=mseq, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                    seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                    pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based)

        return mseq, stats


