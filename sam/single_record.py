"""
Parses Single Sam Record, Sequence Manipulation
"""
import logging
import sam_constants
import align_stats
import Utility
from sam_seq import SamSequence
from sam_constants import SamFlag as SamFlag

LOGGER = logging.getLogger(__name__)

class SamRecord(SamSequence):
    """
    Handles single ended reads or paired reads with missing mate.
    """

    def __init__(self, ref_len, **kwargs):
        """
        Fills in the SamRecord.
        :param int ref_len:  length of reference in nucleotides
        :param dict kwargs: optional keywords.  Valid optional keywords:
        qname, flag, rname, seq, cigar, mapq, qual, pos, rnext, pnext
        :return SamRecord:  an instance of SamRecord
        """
        self.ref_len = ref_len
        self.qname = None
        self.flag = None
        self.rname = None
        self.seq = None
        self.cigar = None
        self.mapq = None
        self.qual = None
        self.pos = None
        self.rnext = None
        self.pnext = None
        self.mate_record = None
        self.nopad_noinsert_seq = None
        self.nopad_noinsert_qual = None
        self.ref_pos_to_insert_seq_qual = None
        self.seq_end_wrt_ref = None
        self.ref_align_len = 0  # includes deletions wrt reference.  Excludes insertions wrt reference
        self.seq_align_len = 0  # includes insertions wrt reference.  Excludes deletions wrt reference

        self.fill_record(**kwargs)


    def is_mapped(self, ref=None):
        """
        :param str ref:  if defined, the record must map to this reference
        :return bool: whether the record is mapped
        """
        if self.flag is None or self.rname is None:
            raise ValueError("No flag, or rname associated with this record")

        return not SamFlag.IS_UNMAPPED & self.flag and (not ref or self.rname == ref)


    def is_mate_mapped(self, ref=None):
        """
        The rname is only set to '*' and pos is only set to '0'
          if all reads in the mate pair are unmapped.
        If the mate is mapped, but this read is not, the rname and pos will be set to the mate's reference hit.

        From SAM specs:
        "For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read
          should have RNAME and POS identical to its mate"

        :param str ref:  if defined, then mate must map to the given reference
        :return bool:  Returns whether the mate for the current record is mapped
        """
        if self.flag is None or self.rnext is None:
            raise ValueError("No flag, rnext associated with this record")

        return (SamFlag.IS_PAIRED & self.flag and not SamFlag.IS_MATE_UNMAPPED & self.flag and
                (not ref or self.rnext == self.rname or self.rnext == "="))


    def is_primary(self):
        """
        :return bool:  Returns whether the current record is primary alignment
        """
        if self.flag is None:
            raise ValueError("No flag associated with this record")
        return not SamFlag.IS_UNMAPPED & self.flag and not SamFlag.IS_SECONDARY_ALIGNMENT & self.flag


    def is_chimeric(self):
        """
        :return bool:  Returns whether the current record is chimeric alignment
        """
        if self.flag is None:
            raise ValueError("No flag associated with this record")
        return SamFlag.IS_CHIMERIC_ALIGNMENT & self.flag


    def is_empty(self):
        """
        Returns whether the sam record has been filled
        :return:
        """
        return not self.qname and not self.seq


    def fill_record(self, **kwargs):
        """
        Fills the sam record using a dict created from a row in a sam file.
        :param dict kwargs:  dict for a row in a sam file.
            Valid keys:  qname, flag, rname, seq, cigar, mapq, qual, pos, rnext, pnext
        """
        if "qname" in kwargs:
            self.qname = kwargs["qname"]
        if "flag" in kwargs:
            self.flag = int(kwargs["flag"])
        if "rname" in kwargs:
            self.rname = kwargs["rname"]
        if "seq" in kwargs:
            self.seq = kwargs["seq"]
        if "cigar" in kwargs:
            self.cigar = kwargs["cigar"]
        if "mapq" in kwargs:
            self.mapq = int(kwargs["mapq"])
        if "qual" in kwargs:
            self.qual = kwargs["qual"]
        if "pos" in kwargs:
            self.pos = int(kwargs["pos"])
        if "rnext" in kwargs:
            self.rnext = kwargs["rnext"]
        if "pnext" in kwargs:
            self.pnext = int(kwargs["pnext"])


        if self.seq == sam_constants.SAM_UNSPECIFIED and self.is_mapped():
            raise ValueError("Sam must specify sequences")

        if self.seq is not None and self.qual is not None and self.qual != sam_constants.SAM_UNSPECIFIED and len(self.seq) != len(self.qual):
            raise ValueError("Expect sequence and quality same length " + str(kwargs) )

        # TODO:  fill in other optional keywords





    def fill_mate(self, mate_record):
        """
        Sets the mate SamRecord.
        :param SamRecord mate_record:  mate's SamRecord
        """
        self.mate_record = mate_record
        if not mate_record.mate_record:
            mate_record.mate_record = self

    def get_name(self):
        """
        :return str:  Return read name as seen in sam file.
        """
        return self.qname


    def get_ref_len(self):
        """
        :return int:  reference length in nucleotides
        """
        return self.ref_len

    def get_slice_intersect_coord(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Gets the 1-based position with respect to the reference of the intersection
        between the sam record sequence and the slice.
        If there is no intersection, returns (None, None)
        :param int slice_start_wrt_ref_1based: 1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based: 1-based slice end position with respect to reference
        :return tuple (int, int): (intersection start, intersection end)
        """
        mate_slice_intersect_start_wrt_ref = None
        mate_slice_intersect_end_wrt_ref = None
        if slice_start_wrt_ref_1based <= self.get_read_end_wrt_ref() and slice_end_wrt_ref_1based >= self.get_read_start_wrt_ref():
            # 1-based position with respect to reference of the start of intersection of the read fragment and slice
            mate_slice_intersect_start_wrt_ref = max(self.get_read_start_wrt_ref(), slice_start_wrt_ref_1based)
            # 1-based position with respect to reference of the end of intersection of read fragment and slice
            mate_slice_intersect_end_wrt_ref = min(self.get_read_end_wrt_ref(), slice_end_wrt_ref_1based)
        return mate_slice_intersect_start_wrt_ref, mate_slice_intersect_end_wrt_ref


    def is_intersect_slice(self, slice_start_wrt_ref_1based, slice_end_wrt_ref_1based):
        """
        Returns whether the read intersects the slice coordinates.
        :param int slice_start_wrt_ref_1based:  1-based slice start position with respect to reference
        :param int slice_end_wrt_ref_1based:  1-based slice end position with respect to reference
        :return:  True if the sequence hits the slice (even with deletions, insertions, or non bases).  False otherwise.
        :rtype: bool
        """
        return slice_start_wrt_ref_1based <= self.get_read_end_wrt_ref() and slice_end_wrt_ref_1based >= self.get_read_start_wrt_ref()


    def get_seq_qual(self, do_pad_wrt_ref=False, do_pad_wrt_slice=False, do_mask_low_qual=False, q_cutoff=10,
                     slice_start_wrt_ref_1based=0, slice_end_wrt_ref_1based=0, do_insert_wrt_ref=False,
                     do_mask_stop_codon=False, stats=None, pad_space_btn_segments="N"):
        """
        Gets the sequence for the sam record.
        If the quality is unspecified in the sam file "*", then returns a string of 0-qual chars the same length as the sequence.
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
        :param str pad_space_btn_segments:  character to use in between aligned segments of a read.  Not used in this method.
        :return tuple (str, str, AlignStats):  (sequence, quality, AlignStats)
        """
        if not stats:
            stats = align_stats.AlignStats()

        # NB:  the original seq from sam file includes softclips.  Remove them with apply_cigar()
        if not self.nopad_noinsert_seq or not self.nopad_noinsert_qual:
            self.__parse_cigar()


        result_seq = ""
        result_qual = ""

        if (not slice_start_wrt_ref_1based and slice_end_wrt_ref_1based) or (slice_start_wrt_ref_1based and not slice_end_wrt_ref_1based):
            raise  ValueError("Either define both slice start and end or don't define either")

        # If not specifed, then the slice is the entire length of the reference
        if not slice_start_wrt_ref_1based:
            slice_start_wrt_ref_1based = 1
        if not slice_end_wrt_ref_1based:
            slice_end_wrt_ref_1based = self.ref_len

        if slice_start_wrt_ref_1based > slice_end_wrt_ref_1based:
            raise ValueError("slice start must be <= slice end")

        read_slice_xsect_start_wrt_ref,  read_slice_xsect_end_wrt_ref= self.get_slice_intersect_coord(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based)
        # Does slice start after the sequence ends or does the slice end before the sequence starts?
        # Then just return empty string or padded gaps wrt ref or slice as desired.
        if not read_slice_xsect_start_wrt_ref and not read_slice_xsect_end_wrt_ref:
            if do_pad_wrt_ref:
                result_seq = SamSequence.do_pad(result_seq, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                          pad_start_wrt_ref=1, pad_end_wrt_ref=self.ref_len, pad_char=sam_constants.SEQ_PAD_CHAR)
                result_qual = SamSequence.do_pad(result_qual, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                          pad_start_wrt_ref=1, pad_end_wrt_ref=self.ref_len, pad_char=sam_constants.QUAL_PAD_CHAR)
            elif do_pad_wrt_slice:
                result_seq = SamSequence.do_pad(result_seq, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                          pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                          pad_char=sam_constants.SEQ_PAD_CHAR)
                result_qual = SamSequence.do_pad(result_qual, seq_start_wrt_ref=None, seq_end_wrt_ref=None,
                                           pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                           pad_char=sam_constants.QUAL_PAD_CHAR)
            return result_seq, result_qual, stats


        # Slice

        # 0-based start position of slice wrt result_seq
        slice_start_wrt_result_seq_0based = read_slice_xsect_start_wrt_ref - self.get_read_start_wrt_ref()
        # 0-based end position of slice wrt result_seq
        slice_end_wrt_result_seq_0based = self.get_ref_align_len() - 1 - (self.get_read_end_wrt_ref() - read_slice_xsect_end_wrt_ref)
        result_seq = self.nopad_noinsert_seq[slice_start_wrt_result_seq_0based:slice_end_wrt_result_seq_0based+1]
        result_qual = self.nopad_noinsert_qual[slice_start_wrt_result_seq_0based:slice_end_wrt_result_seq_0based+1]

        # Mask
        if do_mask_low_qual:
            masked_seq = ""
            masked_qual = ""
            for i, base in enumerate(result_seq):
                base_qual = ord(result_qual[i])-sam_constants.PHRED_SANGER_OFFSET
                if base != sam_constants.SEQ_PAD_CHAR:
                    stats.total_match_1mate += 1
                    if base_qual < q_cutoff:
                        masked_seq += sam_constants.SEQ_PAD_CHAR
                        masked_qual += sam_constants.QUAL_PAD_CHAR
                        stats.total_match_1mate_lo_qual += 1
                    else:
                        masked_seq += base
                        masked_qual += result_qual[i]
                        stats.total_match_1mate_hi_qual += 1
                else:
                    masked_seq += base
                    masked_qual += sam_constants.QUAL_PAD_CHAR

            result_seq = masked_seq
            result_qual = masked_qual

        # Add Insertions
        if do_insert_wrt_ref:
            result_seq_with_inserts = ""
            result_qual_with_inserts = ""
            last_insert_pos_0based_wrt_result_seq = -1  # 0-based position wrt result_seq before the previous insertion

            sliced_insert_dict = self.get_insert_dict(slice_start_wrt_ref_1based, slice_end_wrt_ref_1based)
            # insert_pos_wrt_ref: 1-based reference position before the insertion
            for insert_1based_pos_wrt_ref in sorted(sliced_insert_dict.keys()):  # sort so that insert positions are in order
                insert_seq, insert_qual = sliced_insert_dict[insert_1based_pos_wrt_ref]

                stats.total_insert_blocks += 1
                stats.total_inserts += len(insert_seq)
                stats.total_insert_1mate += len(insert_seq)
                # 0-based position wrt result_seq right before the insertion
                insert_pos_0based_wrt_result_seq =  insert_1based_pos_wrt_ref - read_slice_xsect_start_wrt_ref

                if do_mask_low_qual:
                    masked_insert_seq = ""
                    masked_insert_qual = ""
                    for i, ichar in enumerate(insert_seq):
                        iqual = ord(insert_qual[i])-sam_constants.PHRED_SANGER_OFFSET
                        if iqual >= q_cutoff:  # Only include inserts with high quality
                            masked_insert_seq += ichar
                            masked_insert_qual += insert_qual[i]
                            stats.total_insert_1mate_hi_qual += 1
                        else:
                            stats.total_insert_1mate_lo_qual += 1
                else:
                    masked_insert_seq = insert_seq
                    masked_insert_qual = insert_qual


                result_seq_with_inserts += result_seq[last_insert_pos_0based_wrt_result_seq+1:insert_pos_0based_wrt_result_seq+1] + masked_insert_seq
                result_qual_with_inserts += result_qual[last_insert_pos_0based_wrt_result_seq+1:insert_pos_0based_wrt_result_seq+1] + masked_insert_qual
                last_insert_pos_0based_wrt_result_seq = insert_pos_0based_wrt_result_seq


            result_seq_with_inserts += result_seq[last_insert_pos_0based_wrt_result_seq+1:len(result_seq)]
            result_qual_with_inserts += result_qual[last_insert_pos_0based_wrt_result_seq+1:len(result_qual)]

            if sliced_insert_dict:
                LOGGER.debug("qname=" + self.qname + "insert stats:\n" + stats.dump_insert_stats())

            result_seq = result_seq_with_inserts
            result_qual = result_qual_with_inserts


        # Mask stop codons
        # ASSUME:  that reference starts at beginning of a codon
        if do_mask_stop_codon:
            if read_slice_xsect_start_wrt_ref % Utility.NUC_PER_CODON == 1:
                codon_0based_offset_wrt_result_seq = 0
            else:
                codon_0based_offset_wrt_result_seq = Utility.NUC_PER_CODON - ((read_slice_xsect_start_wrt_ref-1) % Utility.NUC_PER_CODON)

            for nuc_pos_wrt_result_seq_0based in range(codon_0based_offset_wrt_result_seq, len(result_seq), Utility.NUC_PER_CODON):
                codon = result_seq[nuc_pos_wrt_result_seq_0based:nuc_pos_wrt_result_seq_0based+Utility.NUC_PER_CODON]
                if Utility.CODON2AA.get(codon, "") == Utility.STOP_AA:
                    result_seq = result_seq[0:nuc_pos_wrt_result_seq_0based] + "NNN" + result_seq[nuc_pos_wrt_result_seq_0based+Utility.NUC_PER_CODON:]
                    result_qual = result_qual[0:nuc_pos_wrt_result_seq_0based] + (sam_constants.QUAL_PAD_CHAR*3) + result_qual[nuc_pos_wrt_result_seq_0based+Utility.NUC_PER_CODON:]

        # Pad
        if do_pad_wrt_ref:
            result_seq = SamSequence.do_pad(result_seq, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                          seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                          pad_start_wrt_ref=1, pad_end_wrt_ref=self.ref_len, pad_char=sam_constants.SEQ_PAD_CHAR)
            result_qual = SamSequence.do_pad(result_qual, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                           seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                          pad_start_wrt_ref=1, pad_end_wrt_ref=self.ref_len, pad_char=sam_constants.QUAL_PAD_CHAR)

        elif do_pad_wrt_slice:
            result_seq = SamSequence.do_pad(result_seq, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                          seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                          pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                          pad_char=sam_constants.SEQ_PAD_CHAR)
            result_qual = SamSequence.do_pad(result_qual, seq_start_wrt_ref=read_slice_xsect_start_wrt_ref,
                                           seq_end_wrt_ref=read_slice_xsect_end_wrt_ref,
                                           pad_start_wrt_ref=slice_start_wrt_ref_1based, pad_end_wrt_ref=slice_end_wrt_ref_1based,
                                           pad_char=sam_constants.QUAL_PAD_CHAR)



        return result_seq, result_qual, stats


    def get_insert_dict(self, slice_start_wrt_ref_1based=0, slice_end_wrt_ref_1based=0):
        """
        Returns the dict of insert positions to inserts.  Ignores insert positions outside of slice.
        :param int slice_start_wrt_ref_1based: 1-based slice start position with respect to reference.  If 0 or None, uses position 1 wrt ref.
        :param int slice_end_wrt_ref_1based: 1-based slice end position with respect to reference.  If 0 or None, uses reference length.
        :return {int: str}:  dict of 1-based ref position right before the insert => inserted sequence
        """
        if not self.ref_pos_to_insert_seq_qual:
            self.__parse_cigar()
        slice_dict = dict()
        if not slice_start_wrt_ref_1based:
            slice_start_wrt_ref_1based = 1
        if not slice_end_wrt_ref_1based:
            slice_end_wrt_ref_1based = self.get_ref_len()
        for insert_pos_1based, insert in self.ref_pos_to_insert_seq_qual.iteritems():
            if (slice_start_wrt_ref_1based <= insert_pos_1based <= slice_end_wrt_ref_1based-1 or
                    (insert_pos_1based == slice_end_wrt_ref_1based and slice_end_wrt_ref_1based == self.get_ref_len())):
                slice_dict[insert_pos_1based] = insert
        return slice_dict



    def get_read_start_wrt_ref(self):
        """
        Gets the 1-based start position with respect to the reference of the unclipped portion of the sequence.
        :return:
        """
        return self.pos


    def get_read_end_wrt_ref(self):
        """
        Gets the 1-based end position with respect to the reference of the unclipped portion of the sequence.
        :return:
        """
        if not self.seq_end_wrt_ref:
            self.seq_end_wrt_ref =  self.pos + self.get_ref_align_len() - 1
        return self.seq_end_wrt_ref


    def __parse_cigar(self):
        self.nopad_noinsert_seq = ''
        self.nopad_noinsert_qual = ''
        self.ref_pos_to_insert_seq_qual = {}
        self.ref_align_len = 0
        self.seq_align_len = 0
        tokens = sam_constants.CIGAR_RE.findall(self.cigar)

        pos_wrt_seq_0based = 0  # position with respect to sequence, 0based
        pos_wrt_ref_1based = self.pos  # position with respect to reference, 1based
        for token in tokens:
            length = int(token[:-1])
            # Matching sequence: carry it over
            if token[-1] == 'M' or token[-1] == 'X' or token[-1] == '=':
                self.nopad_noinsert_seq += self.seq[pos_wrt_seq_0based:(pos_wrt_seq_0based+length)]
                if self.qual == sam_constants.SAM_UNSPECIFIED:
                    self.nopad_noinsert_qual += sam_constants.QUAL_ZERO_CHAR*length
                else:
                    self.nopad_noinsert_qual += self.qual[pos_wrt_seq_0based:(pos_wrt_seq_0based+length)]
                pos_wrt_seq_0based += length
                pos_wrt_ref_1based += length
                self.ref_align_len += length
                self.seq_align_len += length
            # Deletion relative to reference: pad with gaps
            elif token[-1] == 'D' or token[-1] == 'P' or token[-1] == 'N':
                self.nopad_noinsert_seq += sam_constants.SEQ_PAD_CHAR*length
                self.nopad_noinsert_qual += sam_constants.QUAL_PAD_CHAR*length  # Assign fake placeholder score (Q=-1)
                self.ref_align_len += length
            # Insertion relative to reference: skip it (excise it)
            elif token[-1] == 'I':
                if self.qual == sam_constants.SAM_UNSPECIFIED:
                    # 1-based position wrt reference just BEFORE the insertion
                    self.ref_pos_to_insert_seq_qual.update({pos_wrt_ref_1based - 1:
                                                            (self.seq[pos_wrt_seq_0based:(pos_wrt_seq_0based+length)],
                                                             sam_constants.QUAL_ZERO_CHAR*length)})
                else:
                    # 1-based position wrt reference just BEFORE the insertion
                    self.ref_pos_to_insert_seq_qual.update({pos_wrt_ref_1based - 1:
                                                            (self.seq[pos_wrt_seq_0based:(pos_wrt_seq_0based+length)],
                                                             self.qual[pos_wrt_seq_0based:(pos_wrt_seq_0based+length)])})
                pos_wrt_seq_0based += length
                self.seq_align_len += length
            # Soft clipping leaves the sequence in the SAM - so we should skip it
            elif token[-1] == 'S':
                pos_wrt_seq_0based += length
            elif token[-1] == 'H':  # hard clipping does not leave the sequence in the sam.
                # no-op
                pass
            else:
                raise ValueError("Unable to handle CIGAR token: {} - quitting".format(token))



    def get_ref_align_len(self):
        """
        Returns the length of the alignment with respect to the reference.
        Each inserted base with respect to the reference counts as 0.  Each deleted base with respect to the reference counts as 1.
        :return:
        """
        if not self.ref_align_len:
            self.__parse_cigar()
        return self.ref_align_len


    def get_seq_align_len(self):
        """
        Returns the length of the alignment with respect to the read.
        Each deleted base with respect to the reference counts as 0.  Each inserted base with respect to the reference counts as 1.
        :return:
        """
        if not self.seq_align_len:
            self.__parse_cigar()
        return self.seq_align_len



