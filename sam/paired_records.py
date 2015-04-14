"""
Handles finding merged sequences for paired mates in sam records.
Here, read refers to both mates.
"""

import logging
import sam_constants
import sam.align_stats
import math
from sam.single_record import SamRecord as SamRecord
import Utility
from sam_seq import SamSequence
from collections import OrderedDict

LOGGER = logging.getLogger(__name__)

class PairedRecord(SamSequence):
    """
    Only for mates that hit the same reference.
    """

    def __init__(self, mate1_record, mate2_record):
        self.mate1 = mate1_record
        self.mate2 = mate2_record
        if mate1_record and mate2_record:
            if mate1_record.get_read_start_wrt_ref() > mate2_record.get_read_start_wrt_ref():
                self.mate2 = mate1_record  # switch them so that mate1 is on left, mate2 is on right wrt ref
                self.mate1 = mate2_record

            if self.mate1.rname != self.mate2.rname:
                raise ValueError("mate1_record and mate2_record  must align to the same reference.  " +
                                 "Refs: =(" + str(self.mate1.rname) + ", " + str(self.mate2.rname) + ")")
            if self.mate1.qname != self.mate2.qname:
                raise ValueError("mate1_record and mate2_record should have the same qname.  {}, {}".format(self.mate1.qname, self.mate2.qname))
            if self.mate1.get_ref_len() != self.mate2.get_ref_len():
                raise ValueError("mate1_record and mate2_record should have the same reference length.  {}, {}".format(self.mate1.get_ref_len(), self.mate2.get_ref_len()))
            self.mate1.fill_mate(self.mate2)


        if not self.mate1 and not self.mate2:
            raise ValueError("Must supply at least one mate")
        self.read_start_wrt_ref = None
        self.read_end_wrt_ref = None


    def get_name(self):
        """
        :return str:  Return read name as seen in sam file.
        """
        if self.mate1:
            return self.mate1.get_name()
        else:
            return self.mate2.get_name()

    def get_ref_len(self):
        """
        :return int :
        """
        if self.mate1:
            return self.mate1.get_ref_len()
        else:
            return self.mate2.get_ref_len()

    def get_read_start_wrt_ref(self):
        """
        Gets the 1-based position start of the merged read.
        :return:
        """
        if not self.read_start_wrt_ref:
            if self.mate1 and self.mate2:
                self.read_start_wrt_ref = min(self.mate1.get_read_start_wrt_ref(), self.mate2.get_read_start_wrt_ref())
            elif self.mate1:
                self.read_start_wrt_ref = self.mate1.get_read_start_wrt_ref()
            elif self.mate2:
                self.read_start_wrt_ref = self.mate2.get_read_start_wrt_ref()
        return self.read_start_wrt_ref


    def get_read_end_wrt_ref(self):
        """
        Gets the 1-based position end of the merged read.
        :return:
        """
        if not self.read_end_wrt_ref:
            if self.mate1 and self.mate2:
                self.read_end_wrt_ref = max(self.mate1.get_read_end_wrt_ref(), self.mate2.get_read_end_wrt_ref())
            elif self.mate1:
                self.read_end_wrt_ref = self.mate1.get_read_end_wrt_ref()
            elif self.mate2:
                self.read_end_wrt_ref = self.mate2.get_read_end_wrt_ref()
        return self.read_end_wrt_ref


    def is_intersect_slice(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Returns whether the read intersects the slice coordinates.
        :param int slice_start_wrt_ref_1based:  1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based:  1-based slice end position with respect to reference
        :return:  True if the sequence hits the slice (even with deletions, insertions, or non bases).  False otherwise.
        :rtype: bool
        """
        return slice_start_wrt_ref_1based <= self.get_read_end_wrt_ref() and slice_end_wrt_ref_1based >= self.get_read_start_wrt_ref()


    def is_between_mates(self, pos_wrt_ref_1based):
        """
        Returns whether mates have a gap in between them and the position is inside that gap.
        :param pos_wrt_ref_1based:  1-based position with respect to reference
        :return:
        """
        if not self.mate1 or not self.mate2:
            return False
        return ((self.mate1.get_read_end_wrt_ref() < pos_wrt_ref_1based < self.mate2.get_read_start_wrt_ref()) or
                (self.mate2.get_read_end_wrt_ref() < pos_wrt_ref_1based < self.mate1.get_read_start_wrt_ref()))


    def get_slice_intersect_coord(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Gets the 1-based position with respect to the reference of the intersection
        between the read and the slice.
        If there is no intersection, returns (None, None)
        :param int slice_start_wrt_ref_1based: 1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based: 1-based slice end position with respect to reference
        :return tuple (int, int): (intersection start, intersection end)
        """
        read_slice_intersect_start_wrt_ref = None
        read_slice_intersect_end_wrt_ref = None
        if self.is_intersect_slice(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
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
        if not self.mate1 or not self.mate2:
            return False
        return ((self.mate1.get_read_start_wrt_ref() <= pos_wrt_ref_1based <= self.mate1.get_read_end_wrt_ref()) and
                (self.mate2.get_read_start_wrt_ref() <= pos_wrt_ref_1based <= self.mate2.get_read_end_wrt_ref()))


    @staticmethod
    def calc_q_cutoff_overlap(q_cutoff):
        """
        If overlapping mate bases agree, either mate quality can be lower as long as the
          Probability of single base error at the quality cutoff >= Probability of both mates having same erroneous base
        Quality = -10 log P(error),  P(error) = 1e-(Quality/10)  by Phred Score definition.
        Assume that the probability of a wrong base is independent across mates, positions.
        Say the correct base = G, the permutations of wrong overlapping bases = {AA, AC, AT, CA, CC, CT, TA, TC, TT}.
        ==> There are 3 permutations where the overlapping bases are the same; There are 9 possible permutations of wrong overlapping bases.

        P(single base error at q_cutoff) >= P(both mates same base & both mates wrong base)
        P(single base error at q_cutoff) >= P(both mates wrong base) P(both mates same base | both mates wrong base)
        1e-(q_cutoff/10) >= P(mate1 1-base error) P(mate2 1-base error) [3 matching overlaps / 9 possible permutations of overlapping wrong bases)
        1e-(q_cutoff/10) >= 1e-q1/10 *  1e-q2/10  * 1/3
        3e-(q_cutoff/10) >= 1e-q1/10 *  1e-q2/10
        log3 - q_cutoff/10 >= -q1/10 - q2/10
        -10log3 + qcutoff <= q1 + q2
        qcutoff_overlap = -10log3 + qcutoff <= q1 + q2
        :param int q_cutoff:  single base quality cutoff.
        :return float:  overlapping base quality cutoff, which is much lower than single base quality cutoff.
        """
        q_cutoff_overlap = -10 * math.log10(3) + q_cutoff
        return q_cutoff_overlap



    # TODO:  support other quality formats than those that use Sanger phred offsets
    def __merge_inserts(self, sliced_mseq, sliced_mqual, q_cutoff, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based,
                           stats=None):
        """
        Assumes that sliced_mseq only contains the sliced portion of the merged read sequence.

        :param str sliced_mseq:  merged, sliced read sequence
        :param str sliced_mqual:  merged, sliced quality scores in ASCII.  Assumes Phred Sanger format.
        :param int q_cutoff:  bases with quality lower than this will be masked
        :param int slice_start_wrt_ref_1based: 1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based: 1-based slice end position with respect to reference
        :return str, str, AlignStats:  merged sliced read sequence with inserts, merged sliced qual with inserts, AlignStats
        """
        # track stats about inserts
        if not stats:
            stats = sam.align_stats.AlignStats()

        q_cutoff_overlap = PairedRecord.calc_q_cutoff_overlap(q_cutoff)


        merge_rpos_to_insert = OrderedDict()
        mate1_insert_dict = self.mate1.get_insert_dict(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based) if self.mate1 else dict()
        mate2_insert_dict = self.mate2.get_insert_dict(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based) if self.mate2 else dict()

        # Merge the insert positions from both mates
        uniq_insert_pos = set()
        for x in mate1_insert_dict.keys() + mate2_insert_dict.keys():
            uniq_insert_pos.add(x)

        # For the positions in which the mates overlap, check if the inserts are the same in both mates
        # insert_pos1_wrt_ref_1based is the 1-based position wrt ref before the insertion
        for insert_pos_wrt_ref_1based in sorted(uniq_insert_pos):
            stats.total_insert_blocks += 1
            # Does this insert fall into the mate overlap?
            insert_seq1, insert_qual1 = mate1_insert_dict.get(insert_pos_wrt_ref_1based, ("", ""))
            insert_seq2, insert_qual2 = mate2_insert_dict.get(insert_pos_wrt_ref_1based, ("", ""))
            if self.is_in_mate_overlap(insert_pos_wrt_ref_1based):
                # If the inserted position exists in one mate, then keep it if low quality, else remove
                # If the inserted position exists in both mates, but they disagree on the base, take the higher quality base > q_cutoff.
                # If the inserted position exists in both mates, but they disagree on the base, mask it if both low quality.
                masked_insert_seq = ""
                masked_insert_qual = ""
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
                            masked_insert_qual += (insert_qual1[i] if iqual1 >= iqual2 else insert_qual2[i])
                            stats.total_insert_agree_hi_qual += 1
                        else:
                            masked_insert_seq += "N"
                            masked_insert_qual += sam_constants.QUAL_PAD_CHAR
                            stats.total_insert_agree_lo_qual += 1
                    elif ichar2 and not ichar1 :  # Only mate2 has insert here.  Exclude insert since insert in overlap
                        stats.total_insert_conflict += 1
                        if iqual2 >= q_cutoff:
                            stats.total_insert_conflict_hino_qual += 1
                        else:
                            stats.total_insert_conflict_lono_qual += 1
                    elif ichar1 and not ichar2:  # Only mate1 has insert here.  Exclude insert since insert in overlap
                        stats.total_insert_conflict += 1
                        if iqual1 >= q_cutoff:
                            stats.total_insert_conflict_hino_qual += 1
                        else:
                            stats.total_insert_conflict_lono_qual += 1
                    else:  # Mate1 and Mate2 both have inserts here, but bases conflict.
                        stats.total_insert_conflict += 1

                        if iqual1 < q_cutoff and iqual2 < q_cutoff:  # Both inserts low quality, mask it
                            masked_insert_seq += "N"
                            masked_insert_qual += sam_constants.QUAL_PAD_CHAR
                            stats.total_insert_conflict_lo_qual += 1
                        elif iqual1 > iqual2 >= q_cutoff:  # both inserts high quality, but insert1 higher
                            masked_insert_seq += ichar1
                            masked_insert_qual += insert_qual1[i]
                            stats.total_insert_conflict_hi_qual += 1
                        elif iqual2 > iqual1 >= q_cutoff:  # both sequences high quality, but seq2 higher
                            masked_insert_seq += ichar2
                            masked_insert_qual += insert_qual2[i]
                            stats.total_insert_conflict_hi_qual += 1
                        elif iqual1 >= q_cutoff > iqual2:  # seq1 high quality, seq2 low quality
                            masked_insert_seq += ichar1
                            masked_insert_qual += insert_qual1[i]
                            stats.total_insert_conflict_hilo_qual += 1
                        elif iqual2 >= q_cutoff > iqual1:  # seq2 high quality, seq2 low quality
                            masked_insert_seq += ichar2
                            masked_insert_qual += insert_qual2[i]
                            stats.total_insert_conflict_hilo_qual += 1
                        else:  # both seq high equal quality, mask it:   iqual1 == iqual2 >= q_cutoff
                            masked_insert_seq += "N"
                            masked_insert_qual += sam_constants.QUAL_PAD_CHAR
                            stats.total_insert_conflict_equal_hi_qual += 1


                merge_rpos_to_insert[insert_pos_wrt_ref_1based] = (masked_insert_seq, masked_insert_qual)
            else:
                # If the the insert is in a region that does not overlap with the mate, exclude it if it's low quality
                masked_insert_seq = ""
                masked_insert_qual = ""
                insert_seq = insert_seq1 if insert_seq1 else insert_seq2
                insert_qual = insert_qual1 if insert_qual1 else insert_qual2
                for i, ichar in enumerate(insert_seq):
                    stats.total_insert_1mate += 1
                    stats.total_inserts += 1
                    iqual = ord(insert_qual[i])-sam_constants.PHRED_SANGER_OFFSET
                    if iqual >= q_cutoff:
                        masked_insert_seq += ichar
                        masked_insert_qual += insert_qual[i]
                        stats.total_insert_1mate_hi_qual += 1
                    else:
                        stats.total_insert_1mate_lo_qual += 1

                merge_rpos_to_insert[insert_pos_wrt_ref_1based] = (masked_insert_seq, masked_insert_qual)




        mseq_with_inserts = ""
        mqual_with_inserts = ""
        last_insert_pos_0based_wrt_mseq = -1  # 0-based position wrt result_seq before the previous insertion
        # insert_pos_wrt_ref: 1-based reference position before the insertion
        for insert_1based_pos_wrt_ref, (insert_seq, insert_qual) in merge_rpos_to_insert.iteritems():
            # 0-based position wrt sliced_mseq right before the insertion
            insert_pos_0based_wrt_mseq =  insert_1based_pos_wrt_ref - slice_start_wrt_ref_1based
            mseq_with_inserts += sliced_mseq[last_insert_pos_0based_wrt_mseq+1:insert_pos_0based_wrt_mseq+1] + insert_seq
            mqual_with_inserts += sliced_mqual[last_insert_pos_0based_wrt_mseq+1:insert_pos_0based_wrt_mseq+1] + insert_qual

            last_insert_pos_0based_wrt_mseq = insert_pos_0based_wrt_mseq

        mseq_with_inserts += sliced_mseq[last_insert_pos_0based_wrt_mseq+1:len(sliced_mseq)]
        mqual_with_inserts += sliced_mqual[last_insert_pos_0based_wrt_mseq+1:len(sliced_mqual)]

        if merge_rpos_to_insert:
            LOGGER.debug("qname=" + self.get_name() + " insert stats:\n" + stats.dump_insert_stats())

        return mseq_with_inserts, mqual_with_inserts, stats



    def __merge_match(self, q_cutoff=10, pad_space_btn_mate="N", slice_start_wrt_ref_1based=0, slice_end_wrt_ref_1based=0, stats=None):
        """
        Merge bases in reads that align as a match to the reference (as opposed to bases that align as an insert to the ref).
        Does not do validation of slice coordinates.
        :param int q_cutoff:  bases with quality lower than this will be masked
        :param int slice_start_wrt_ref_1based: 1-based slice start position with respect to reference.
        :param int slice_end_wrt_ref_1based: 1-based slice end position with respect to reference.
        :param AlignStats stats:  Keeps track of AlignStats.  Creates a new instance if None.
        :return str, str:  merged sequence, merged quality.  When there is discordant bases, the quality score for the winning base is chosen.
        The quality score for deleted or masked bases is set to zero.
        """
        mseq = ""
        mqual = ""

        if not stats:
            stats = sam.align_stats.AlignStats()

        # Extract portion of reads that fit within the intersection of the read fragment & slice.
        # Pad the extracted sequences just enough so that they line up within the intersection.
        # Do not pad with respect to the reference yet.  Do not pad with respect to the slice yet.
        # Pop out insertions with respect to the reference so that it is easier to find the sequence coordinates wrt ref.
        # Keep track of the insertions and the ref pos right before the insertion.
        seq1, qual1 = "", ""
        seq2, qual2 = "", ""
        if self.mate1:
            seq1, qual1, stats = self.mate1.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True,
                                                         do_mask_low_qual=False, q_cutoff=q_cutoff,
                                                         slice_start_wrt_ref_1based=slice_start_wrt_ref_1based,
                                                         slice_end_wrt_ref_1based=slice_end_wrt_ref_1based,
                                                         do_insert_wrt_ref=False, do_mask_stop_codon=False, stats=stats)

            if not self.mate2: # in case one mate is defined but not the other
                seq2 = sam_constants.SEQ_PAD_CHAR * len(seq1)
                qual2 = sam_constants.QUAL_PAD_CHAR * len(qual1)


        if self.mate2:
            seq2, qual2, stats = self.mate2.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=True,
                                                     do_mask_low_qual=False, q_cutoff=q_cutoff,
                                                     slice_start_wrt_ref_1based=slice_start_wrt_ref_1based,
                                                     slice_end_wrt_ref_1based=slice_end_wrt_ref_1based,
                                                     do_insert_wrt_ref=False, do_mask_stop_codon=False, stats=stats)

            if not self.mate1:  # in case one mate is defined but not the other
                seq1 = sam_constants.SEQ_PAD_CHAR * len(seq2)
                qual1 = sam_constants.QUAL_PAD_CHAR * len(qual2)


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
                mqual += sam_constants.QUAL_PAD_CHAR
            # only mate2 has real base here
            elif seq1[i] == sam_constants.SEQ_PAD_CHAR and seq2[i] != sam_constants.SEQ_PAD_CHAR:
                stats.total_match_1mate += 1
                if q2 >= q_cutoff:
                    mseq += seq2[i]
                    mqual += qual2[i]
                    stats.total_match_1mate_hi_qual += 1
                else:
                    mseq += "N"
                    mqual += sam_constants.QUAL_PAD_CHAR
                    stats.total_match_1mate_lo_qual += 1

            # only mate1 has real base here
            elif seq2[i] == sam_constants.SEQ_PAD_CHAR and seq1[i] != sam_constants.SEQ_PAD_CHAR:
                stats.total_match_1mate += 1
                if q1 >= q_cutoff:
                    mseq += seq1[i]
                    mqual += qual1[i]
                    stats.total_match_1mate_hi_qual += 1
                else:
                    mseq += "N"
                    mqual += sam_constants.QUAL_PAD_CHAR
                    stats.total_match_1mate_lo_qual += 1

            # Both mates have the same base at this position
            elif seq1[i] == seq2[i] and seq1[i] != sam_constants.SEQ_PAD_CHAR:
                stats.total_match_nonconflict += 1
                if q_cutoff_overlap <= q1 + q2:
                    mseq += seq1[i]
                    mqual += (qual1[i] if q1 >= q2 else qual2[i])
                    stats.total_match_nonconflict_hi_qual += 1
                else:
                    mseq += "N"
                    mqual += sam_constants.QUAL_PAD_CHAR
                    stats.total_match_nonconflict_lo_qual += 1

            # Both mates have disagreeing bases at this position: take the high confidence
            elif seq1[i] != seq2[i]:
                stats.total_match_conflict += 1
                if q1 < q_cutoff and q2 < q_cutoff:  # both sequences low quality
                    mseq += "N"
                    mqual += sam_constants.QUAL_PAD_CHAR
                    stats.total_match_conflict_lo_qual += 1
                elif q1 > q2 >= q_cutoff:  # both sequences high quality, but seq1 higher
                    mseq += seq1[i]
                    mqual += qual1[i]
                    stats.total_match_conflict_hi_qual += 1
                elif q1 >= q_cutoff > q2:  # seq1 high quality, seq2 low quality
                    mseq += seq1[i]
                    mqual += qual1[i]
                    stats.total_match_conflict_hilo_qual += 1
                elif q2 > q1 >= q_cutoff:  # both sequences high quality, but seq2 higher
                    mseq += seq2[i]
                    mqual += qual2[i]
                    stats.total_match_conflict_hi_qual += 1
                elif q2 >= q_cutoff > q1:  # seq2 high quality, seq2 low quality
                    mseq += seq2[i]
                    mqual += qual2[i]
                    stats.total_match_conflict_hilo_qual += 1
                elif q1 == q2 >= q_cutoff:
                    mseq += "N"
                    mqual += sam_constants.QUAL_PAD_CHAR
                    stats.total_match_conflict_equal_hi_qual += 1
                else:
                    raise ValueError("We should never get here.  Unanticipated use case for merging sam reads")

            else:
                raise ValueError("We should never get here.  Unanticipated use case for merging sam reads")

        return mseq, mqual


    def get_seq_qual(self, q_cutoff=10, pad_space_btn_segments="N", do_insert_wrt_ref=False, do_pad_wrt_ref=True,
                        do_pad_wrt_slice=False, do_mask_stop_codon=False, slice_end_wrt_ref_1based=0, slice_start_wrt_ref_1based=0,
                        stats=None):
        """
        Merge two sequences that overlap over some portion (paired-end
        reads).  Using the positional information in the SAM file, we will
        know where the sequences lie relative to one another.  In the case
        that the base in one read has no complement in the other read
        (in partial overlap region), take that base at face value.

        Allows insertions.


        When merging quality, the quality score for the winning base is chosen for discordant bases.
        The quality score for deleted or masked bases is set to zero.

        :param q_cutoff:
        :param do_insert_wrt_ref:
        :param do_pad_wrt_ref:
        :param do_mask_stop_codon:
        :param slice_end_wrt_ref_1based:
        :param slice_start_wrt_ref_1based:
        :param stats:
        :param pad_space_btn_segments:
        :param do_pad_wrt_slice:
        :param int q_cutoff: quality cutoff below which a base is converted to N if there is no consensus between the mates.
        :param bool do_insert_wrt_ref: whether insertions with respect to reference is allowed.
        :param bool do_mask_stop_codon: If True, then masks stop codons with "NNN".  Assumes that reference starts at beginning of codon.
                Performs stop codon masking after masking low quality bases and after including insertions (if do_insert_wrt_ref==True).
        :return:  merged paired-end sequence, merged quality, AlignStats instance
        :rtype : (str, str, AlignStats)
        """

        mseq = ""
        mqual = ""
        if not stats:
            stats = sam.align_stats.AlignStats()

        if not slice_start_wrt_ref_1based:
            slice_start_wrt_ref_1based = 1
        if not slice_end_wrt_ref_1based:
            slice_end_wrt_ref_1based = self.get_ref_len()
        if slice_start_wrt_ref_1based > slice_end_wrt_ref_1based:
            raise ValueError("slice start must be before slice end")

        # If the slice does not intersect either mate,
        # Then just return empty string or padded gaps wrt ref or slice as desired.
        if not self.is_intersect_slice(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
            if do_pad_wrt_ref:
                mseq = SamSequence.do_pad(mseq, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                        pad_start_wrt_ref=1, pad_end_wrt_ref=self.mate1.ref_len,
                                        pad_char=sam_constants.SEQ_PAD_CHAR)
                mqual = SamSequence.do_pad(mqual, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                        pad_start_wrt_ref=1, pad_end_wrt_ref=self.mate1.ref_len,
                                        pad_char=sam_constants.QUAL_PAD_CHAR)
            elif do_pad_wrt_slice:
                mseq = SamSequence.do_pad(mseq, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                        pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                        pad_char=sam_constants.SEQ_PAD_CHAR)
                mqual = SamSequence.do_pad(mqual, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                        pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                        pad_char=sam_constants.QUAL_PAD_CHAR)
            return mseq, mqual, stats


        # 1-based position with respect to reference of intersection between the read fragment and slice
        read_slice_xsect_start_wrt_ref, read_slice_xsect_end_wrt_ref = self.get_slice_intersect_coord(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based)

        mseq, mqual = self.__merge_match(q_cutoff=q_cutoff, pad_space_btn_mate=pad_space_btn_segments,
                         slice_start_wrt_ref_1based=read_slice_xsect_start_wrt_ref,
                         slice_end_wrt_ref_1based=read_slice_xsect_end_wrt_ref, stats=stats)


        if do_insert_wrt_ref:
            mseq, mqual, stats = self.__merge_inserts(sliced_mseq=mseq, sliced_mqual=mqual, q_cutoff=q_cutoff,
                                                  slice_start_wrt_ref_1based=read_slice_xsect_start_wrt_ref,
                                                  slice_end_wrt_ref_1based=read_slice_xsect_end_wrt_ref, stats=stats)

        # Mask stop codons
        # ASSUME:  that reference starts at beginning of a codon
        # TODO:  handle situation if inserts before slice cause extra or remove stop codon
        if do_mask_stop_codon:
            if read_slice_xsect_start_wrt_ref % Utility.NUC_PER_CODON == 1:
                codon_0based_offset_wrt_result_seq = 0
            else:
                codon_0based_offset_wrt_result_seq = Utility.NUC_PER_CODON - ((read_slice_xsect_start_wrt_ref-1) % Utility.NUC_PER_CODON)

            for nuc_pos_wrt_result_seq_0based in range(codon_0based_offset_wrt_result_seq, len(mseq), Utility.NUC_PER_CODON):
                codon = mseq[nuc_pos_wrt_result_seq_0based:nuc_pos_wrt_result_seq_0based+Utility.NUC_PER_CODON]
                if Utility.CODON2AA.get(codon, "") == Utility.STOP_AA:
                    mseq = mseq[0:nuc_pos_wrt_result_seq_0based] + "NNN" + mseq[nuc_pos_wrt_result_seq_0based+Utility.NUC_PER_CODON:]
                    mqual = mqual[0:nuc_pos_wrt_result_seq_0based] + (sam_constants.QUAL_PAD_CHAR*3) + mqual[nuc_pos_wrt_result_seq_0based+Utility.NUC_PER_CODON:]

        # Now pad with respect to reference
        if do_pad_wrt_ref:
            mseq = SamSequence.do_pad(seq=mseq, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                    seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                    pad_start_wrt_ref=1, pad_end_wrt_ref=self.mate1.ref_len, pad_char=sam_constants.SEQ_PAD_CHAR)
            mqual = SamSequence.do_pad(seq=mqual, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                    seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                    pad_start_wrt_ref=1, pad_end_wrt_ref=self.mate1.ref_len, pad_char=sam_constants.QUAL_PAD_CHAR)

        elif do_pad_wrt_slice:
            mseq = SamSequence.do_pad(seq=mseq, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                    seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                    pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                    pad_char=sam_constants.SEQ_PAD_CHAR)
            mqual = SamSequence.do_pad(seq=mqual, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                    seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                    pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                    pad_char=sam_constants.QUAL_PAD_CHAR)
        return mseq, mqual, stats


