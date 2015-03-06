#  Keeps track of information for testing sam parsing
import re


class SamTestCase:

    class __MateInfo:
        """
        Helper class for holding mate info
        """
        def __init__(self, qname, flag, rname, pos, mapq, cigar, tlen, seq, qual):
            self.qname = qname
            self.flag = flag
            self.rname = rname
            self.pos = pos
            self.mapq = mapq
            self.cigar = cigar
            self.rnext = "*"
            self.pnext = 0
            self.tlen = tlen
            self.seq = seq
            self.qual = qual




    """
    Represents a testcase for merging paired reads.  Only handles situtations in which each mate has 0 or 1 alignment.
    Doesn't handle chimeric or secondary alignments.
    """
    def __init__(self, test_desc, target_ref, ref2len, slice_start, slice_end):
        self.test_desc = test_desc
        self.target_ref = target_ref
        self.ref2len = ref2len
        self.slice_start = slice_start
        self.slice_end = slice_end
        self.read_name = ""     # Should be the same for both mates
        self.merged_read_seq = ""
        self.mate1 = None
        self.mate2 = None
        self.insert_pos_wrt_merged_seq_1based = []




    def add_read(self, is_first, qname, flag, rname, mapq, cigar, tlen, seq, qual):
        """
        Adds a mate into the testcase.  Expects that the mates are multiple-sequence aligned to each other.
        Stores the mate sequences and qualities in unaligned form.
        :param is_first:
        :param qname:
        :param flag:
        :param rname:
        :param mapq:
        :param cigar:
        :param tlen:
        :param seq:
        :param qual:
        :return:
        """
        if len(seq) != len(qual):
            raise ValueError("seq should be same length as qual: " + seq + " " + qual)

        unpad_seq = seq.replace("-", "")
        unpad_qual = qual.replace("-", "")
        if len(unpad_seq) != len(unpad_qual):
            raise ValueError("ungapped seq should be same length as ungapped qual: " + + seq + " " + qual)
        self.read_name = qname
        first_non_gap = re.search(r"[^-]", seq)
        pos = first_non_gap.start() + 1
        if is_first:
            if self.mate2  and self.mate2.qname != qname:
                raise ValueError("Mates must have same name, mate1=" + qname  + " mate2=" + self.mate2.qname)
            self.mate1 = SamTestCase.__MateInfo(qname, flag, rname, pos, mapq, cigar, tlen, unpad_seq, unpad_qual)
            if self.mate2:
                self.mate1.rnext = self.mate2.rname
                self.mate1.pnext = self.mate2.pos
        else:
            if self.mate1  and self.mate1.qname != qname:
                raise ValueError("Mates must have same name, mate2=" + qname  + " mate1=" + self.mate1.qname)
            self.mate2 = SamTestCase.__MateInfo(qname, flag, rname, pos, mapq, cigar, tlen, unpad_seq, unpad_qual)
            if self.mate1:
                self.mate2.rnext = self.mate1.rname
                self.mate2.pnext = self.mate1.pos


    def add_merged_read(self, seq, inserts):
        self.merged_read_seq = seq
        self.insert_pos_wrt_merged_seq_1based = inserts


    def is_same_ref(self):
        """
        :return: whether the mates hit the same ref for their primary alignment
        :rtype:  bool
        """
        return self.mate1 and self.mate2 and self.mate1.rname == self.mate2.rname


    def create_sam_lines(self):
        """
        Creates sam alignment records for the testcase
        :return str  : string representation of newline separated sam records
        """
        lines = ""
        for mate in [self.mate1, self.mate2]:
            if mate:
                lines += "\t".join([mate.qname,
                                    str(mate.flag),
                                    mate.rname,
                                    str(mate.pos),
                                    str(mate.mapq),
                                    mate.cigar,
                                    mate.rnext,
                                    str(mate.pnext),
                                    str(mate.tlen),
                                    mate.seq,
                                    mate.qual]) + "\n"
        return lines


    def get_sliced_merged_read(self, slice_start_pos_wrt_ref_1based=None, slice_end_pos_wrt_ref_1based=None,
                               do_pad_wrt_slice=True, do_insert_wrt_ref=False, do_mask_stop_codon=False):
        """
        Slices the merged sequence at the coordinates wrt ref.
        ASSUMES:
        Merged sequence is already padded wrt ref.


        :param int slice_start_pos_wrt_ref_1based:  1-based start position with respect to reference.  If None, then uses full read.
        :param int slice_end_pos_wrt_ref_1based:  1-based end position with respect to reference.  If None, then uses full read.
        :param bool do_pad_wrt_slice:  If False, then strips left and right pad gaps.
        :param bool do_insert_wrt_ref:  If False, then strips inserts with respect to reference.
        :return: Merged sequence of the paired mates, sliced at desired region.
        :rtype:  str
        """
        slice_start_pos_wrt_ref_1based = 1 if not slice_start_pos_wrt_ref_1based else slice_start_pos_wrt_ref_1based
        slice_end_pos_wrt_ref_1based = self.ref2len[self.target_ref] if not slice_end_pos_wrt_ref_1based else slice_end_pos_wrt_ref_1based

        # Start with the full merged read sequence that includes insertions wrt ref, stop codons, left and right padding wrt reference
        # Strip out inserts if necessary
        stripped_insert_seq = self.merged_read_seq
        if not do_insert_wrt_ref:
            stripped_insert_seq = ""
            for pos_wrt_merged_seq_0based, seqchar in enumerate(self.merged_read_seq):
                pos_wrt_merged_seq_1based = pos_wrt_merged_seq_0based + 1
                if pos_wrt_merged_seq_1based not in self.insert_pos_wrt_merged_seq_1based:
                    stripped_insert_seq += seqchar


        # Mask stop codons
        stop_codon_masked_seq = stripped_insert_seq
        if do_mask_stop_codon:
            for pos in range(0, len(stripped_insert_seq), 3):
                codon = stripped_insert_seq[pos:pos+3]
                if codon in ["TAA", "TGA", "TAG"]:
                    stop_codon_masked_seq =  stripped_insert_seq[0: pos] + "NNN" + stripped_insert_seq[pos+3:]


        # Slice
        sliced_merged_seq = stop_codon_masked_seq[slice_start_pos_wrt_ref_1based-1:slice_end_pos_wrt_ref_1based]
        if do_insert_wrt_ref:  #If there are inserts wrt ref, then we need to translate positions wrt to seq to positions wrt ref
            pos_wrt_ref = 0
            sliced_merged_seq = ""
            for pos_wrt_merged_seq_0based, seqchar in enumerate(stop_codon_masked_seq):
                pos_wrt_merged_seq_1based = pos_wrt_merged_seq_0based + 1
                if pos_wrt_merged_seq_1based not in self.insert_pos_wrt_merged_seq_1based:
                    pos_wrt_ref += 1
                if slice_start_pos_wrt_ref_1based <= pos_wrt_ref <= slice_end_pos_wrt_ref_1based:
                    sliced_merged_seq += seqchar
                elif pos_wrt_ref > slice_end_pos_wrt_ref_1based:
                    break

        if not do_pad_wrt_slice:
            sliced_merged_seq = sliced_merged_seq.lstrip("-").rstrip("-")


        return sliced_merged_seq


    @staticmethod
    def parse_fastq(fastq):
        """
        Parses a fastq file for testcases for merging, masking low quality bases and conflicts, & padding sequences from SAM records.

        Fastq format is used instead of direct SAM or some other format since it's easier for the tester to
            visualize sequence alignments in fastq format.

        There can be a single target reference in the fastq, but multiple references can be listed.

        Each testcase is represented by a testcase description comment, exactly 2 mates of a read,  and the expected merged sequence.
        The expected merged sequence should have low quality bases masked, conflicts masked or excluded as appropriate,
            padded with respect to the reference.  It will include all insertions but will not mask stop codons.
            In the header, there should be an inserts=<> field in the header
            with a comma-delimited list of 1-based positions in the merged sequence of any insertions with respect to the reference.
            Use SamTestCase.get_sliced_merged() to slice the expected merged sequence, strip out inserts, and mask stop codons.

        All sequences and qualities are multiple sequence aligned.
        Each read should include the sam fields flag, rname, mapq, cigar, tlen in the header.
        The pos, rnext, pnext will be inferred by parse_fastq based on the multiple sequence alignment with the reference and the cigar string.

        EG)  There is 1 target reference, "ref1" and 1 testcase
        @ref1   Target  Slice=4,13
        ACGTACGTGGGG
        +
        HHHHHHHHHHHH
        # Comment describing testcase.
        @read<testcase num>/1   flag=99,rname=ref1,mapq=40,cigar=8M,tlen=10
        ACGTACGT----
        +
        HHHHHHHH-----
        @read<testcase num>/2   flag=147,rname=ref1,mapq=40,cigar=8M,tlen=10
        --------GGGG
        +
        --------HHHH
        @merged_read<testcase num>  inserts=<comma separated 1-based positions with respect to merged sequence of inserts with respect to the reference>
        ACGTACGTGGGG
        +
        HHHHHHHHHHHH


        ASSUMES:
        - ORDER:  Reference Sequences are listed first.  Then reads.  For each read, the order is 1st mate, 2nd mate, then expected merged sequence for both mates.
        - References are listed first.  Reference Header Format:  @ref<number>
        - If the reference is the target reference, then it should have header format:  @ref<number>\tTarget\tSlice=<slice start>,<slice_end>
            Merged sequences will be tested for slicing only at the slice defined in the target reference header.
        - Read header format:  read<number>      rname=<ref>,cigar=<cigar>,rnext=<rnext>,tlen=<tlen>
        - Merged read header format:  merged_read<number> <comments about how merged read generated>
        :return :  list of MergeTestCase
        :rtype:  [MergeTestCase]
        """
        merge_tests = []

        target_ref = None
        ref2len = {"*":0}  # reference name : reference sequence length
        slice_start = None
        slice_end = None
        qname = None
        with open(fastq, 'rU') as fh_in:
            merge_test = None
            for i, line in enumerate(fh_in):
                line = line.rstrip()
                if line.startswith("@ref"):
                    line_arr = line[1:].split()
                    ref_name = line_arr[0]
                    ref_seq = fh_in.next().rstrip()
                    sep = fh_in.next()
                    ref_qual = fh_in.next()  # we don't care about the quality of the reference
                    ref2len[ref_name] = len(ref_seq)

                    if len(line_arr) == 3 and line_arr[1] == "Target" and line_arr[2].startswith("Slice="):
                        target_ref = ref_name
                        slice_start, slice_end = [int(x) for x in line_arr[2].replace("Slice=", "").split(",")]

                elif line.startswith("#"):
                    #1 line for test case description
                    test_desc = line[1:]
                    merge_test = SamTestCase(test_desc, target_ref, ref2len, slice_start, slice_end)
                elif line.startswith("@read"):
                    # format:
                    # @read1/1    flag=99,rname=ref1,mapq=40,cigar=4M1I4M,rnext=ref1,tlen=10
                    qname, comments = line.split()
                    is_first =  qname.endswith("/1")
                    qname = qname[1:-2]  # remove @.  remove /1
                    flag, rname, mapq, cigar, tlen = comments.split(",")
                    flag = flag.replace("flag=","")
                    rname = rname.replace("rname=","")
                    mapq = mapq.replace("mapq=","")
                    cigar = cigar.replace("cigar=","")
                    tlen = tlen.replace("tlen=","")

                    seq = fh_in.next().rstrip()
                    sep = fh_in.next()
                    qual = fh_in.next().rstrip()

                    if len(seq) != len(qual):
                        raise ValueError("seq should be same length as qual: line :" + line)
                    if len(seq.lstrip("-").rstrip("-")) != len(qual.lstrip("-").rstrip("-")):
                        raise ValueError("ungapped seq should be same length as ungapped qual:  line :" + line)

                    merge_test.add_read(is_first=is_first, qname=qname, flag=flag, rname=rname, mapq=mapq, cigar=cigar,
                                        tlen=tlen, seq=seq, qual=qual)

                elif line.startswith("@merged_read"):
                    if len(line[1:].split()) > 1:
                        merged_read_name, insert_str = line[1:].split()
                        inserts = [int(x) for x in insert_str.replace("inserts=", "").split(",")]
                    else:
                        merged_read_name = line[1:]
                        inserts = []

                    if merged_read_name != "merged_" + qname:
                        raise ValueError("Something wrong with the read names:  line :" + line)

                    merged_read_seq = fh_in.next().rstrip()
                    sep = fh_in.next()
                    merged_read_qual = fh_in.next().rstrip()
                    i += 3
                    merge_test.add_merged_read(merged_read_seq, inserts)
                    merge_tests.append(merge_test)
                else:
                    raise ValueError("Shouldn't get here, invalid test fastq format. line " + str(i+1) + ":" + line)

        return merge_tests
