import unittest
import sam.paired_records as pair
from sam.sam_record import SamRecord
import re
import csv
import sam.sam_constants as sam_constants
import logging
import sys
import sam.align_stats

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

class MergeTestCase:
    def __init__(self, test_desc):
        self.test_desc = test_desc
        self.merged_read_name = ""
        self.merged_read_seq = ""
        self.mate1 = dict()
        self.mate2 = dict()
        self.insert_pos_wrt_merged_seq_1based = []


    def add_read(self, is_first, qname, flag, rname, mapq, cigar, rnext, tlen, seq, qual):
        self.merged_read_name = qname
        first_non_gap = re.search(r"[^-]", seq)
        pos = first_non_gap.start() + 1
        read_dict = {"qname": qname,
                     "flag": flag,
                     "rname": rname,
                     "pos": pos,
                     "mapq": mapq,
                     "cigar": cigar,
                     "rnext": rnext,
                     "tlen": tlen,
                     "seq": seq,
                     "qual": qual}
        if is_first:
            if self.mate2.get("qname")  and self.mate2["qname"] != qname:
                raise ValueError("Mates must have same name, mate1=" + qname  + " mate2=" + self.mate2["qname"])
            self.mate1 = read_dict
            self.mate1["pnext"] = self.mate2.get("pos", 0)
            self.mate2["pnext"] = pos

        else:
            if self.mate1.get("qname")  and self.mate1["qname"] != qname:
                raise ValueError("Mates must have same name, mate2=" + qname  + " mate1=" + self.mate1["qname"])
            self.mate2 = read_dict
            self.mate2["pnext"] = self.mate1.get("pos", 0)
            self.mate1["pnext"] = pos


    def add_merged_read(self, seq, inserts):
        self.merged_read_seq = seq
        self.insert_pos_wrt_merged_seq_1based = inserts


    def create_sam_lines(self, is_include_unalign=False):
        """
        Creates sam alignment entries for the reads
        :param bool    is_include_unalign:  whether to include records for unaligned mates
        :return [str]  sam alignment lines:
        """
        lines = ""
        for mate in [self.mate1, self.mate2]:
            if is_include_unalign or mate:
                if not int(mate["flag"]) & sam_constants.SamFlag.IS_UNMAPPED and mate["rname"] == "ref1":
                    lines += "\t".join([mate["qname"],
                                       str(mate["flag"]),
                                       mate["rname"],
                                       str(mate["pos"]),
                                       str(mate["mapq"]),
                                       mate["cigar"],
                                       mate["rnext"],
                                       str(mate["pnext"]),
                                       str(mate["tlen"]),
                                       # Remove left and right padding
                                       # it's only for our visual benefit in the fastq and shouldn't be in sam
                                       mate["seq"].lstrip("-").rstrip("-"),
                                       mate["qual"].lstrip("-").rstrip("-")
                                       ]) + "\n"
        return lines


    def get_sliced_merged_read(self, slice_start_pos_wrt_ref_1based, slice_end_pos_wrt_ref_1based):
        """
        Slices the merged sequence at the coordinates wrt ref.
        ASSUMES:
        Merged sequence is always padded wrt ref.
        :param slice_start_pos_wrt_ref_1based:
        :param slice_end_pos_wrt_ref_1based:
        :return:
        """
        pos_wrt_ref = 0
        slice_start_pos_wrt_merged_seq_0based = None
        slice_end_pos_wrt_merged_seq_0based = None
        for pos_wrt_merged_seq_0based, seqchar in enumerate(self.merged_read_seq):
            pos_wrt_merged_seq_1based = pos_wrt_merged_seq_0based + 1
            if pos_wrt_merged_seq_1based not in self.insert_pos_wrt_merged_seq_1based:
                pos_wrt_ref += 1
                if pos_wrt_ref == slice_start_pos_wrt_ref_1based:
                    slice_start_pos_wrt_merged_seq_0based = pos_wrt_merged_seq_0based
                if pos_wrt_ref == slice_end_pos_wrt_ref_1based:
                    slice_end_pos_wrt_merged_seq_0based  = pos_wrt_merged_seq_0based


        if not slice_start_pos_wrt_merged_seq_0based or not slice_start_pos_wrt_merged_seq_0based:
            raise  ValueError("Invalid slice calculations")
        return self.merged_read_seq[slice_start_pos_wrt_merged_seq_0based:slice_end_pos_wrt_merged_seq_0based+1]






class TestSliceAligned(unittest.TestCase):

    #TODO:  testcases:
    # - merge reads where 1 mate aligns to a different reference
    # - merge reads where 1 mate dpesn't align at all

    REF_LEN = 0
    TEST_PAD_INSERT_FASTQ = "./data/test.pad.insert.fq"  # Test reads for allowing inserts and padding wrt reference
    TEST_SLICE_INSERT_FASTQ = "./data/test.slice.insert.fq"  # Test reads for allowing inserts and padding wrt reference
    Q_CUTOFF = 20
    M_CUTOFF = 20

    def __generate_merged_read_sam(self, fastq):
        """
        From a fastq file, generates the sam file for testing paired-end sequence merging.
        Comments are allowed and describe the following testcase.

        ASSUMES:
        - ORDER:  Reference Sequences are listed first.  Then reads.  For each read, the merged read follows.
        - References are listed first.  Format:  ref<number>
        - Read header format:  read<number>      rname=<ref>,cigar=<cigar>,rnext=<rnext>,tlen=<tlen>
        - Merged read header format:  merged_read<number> <comments about how merged read generated>
        :return :  tuple of sam filename and list of MergeTestCase
        :rtype:  (str, [MergeTestCase])
        """
        merge_tests = []
        sam = fastq.replace(".fq", ".sam")

        REF_LEN = None
        REFERENCE = None
        with open(fastq, 'rU') as fh_in, open(sam, 'w') as fh_out:
            fh_out.write("@HD\tVN:1.4\tSO:queryname]\n")
            fh_out.write("@SQ\tSN:consensus\tLN:10\n")

            merge_test = None
            for i, line in enumerate(fh_in):
                line = line.rstrip()
                if line.startswith("@ref"):
                    ref_name = line[1:]
                    ref_seq = fh_in.next().rstrip()
                    sep = fh_in.next()
                    ref_qual = fh_in.next().rstrip()
                    REF_LEN = len(ref_seq)
                    REFERENCE = ref_name
                elif line.startswith("#"):
                    #1 line for test case description
                    test_desc = line[1:]
                    merge_test = MergeTestCase(test_desc)
                elif line.startswith("@read"):
                    # format:
                    # @read1/1    flag=99,rname=ref1,mapq=40,cigar=4M1I4M,rnext=ref1,tlen=10
                    qname, comments = line.split()
                    is_first =  qname.endswith("/1")
                    qname = qname[1:-2]  # remove @.  remove /1
                    flag, rname, mapq, cigar, rnext, tlen = comments.split(",")
                    flag = flag.replace("flag=","")
                    rname = rname.replace("rname=","")
                    mapq = mapq.replace("mapq=","")
                    cigar = cigar.replace("cigar=","")
                    rnext = rnext.replace("rnext=","")
                    tlen = tlen.replace("tlen=","")

                    seq = fh_in.next().rstrip()
                    sep = fh_in.next()
                    qual = fh_in.next().rstrip()
                    if len(seq) != len(qual):
                        raise ValueError("seq should be same length as qual: " + line)
                    if len(seq.lstrip("-").rstrip("-")) != len(qual.lstrip("-").rstrip("-")):
                        raise ValueError("ungapped seq should be same length as ungapped qual: " + line)

                    merge_test.add_read(is_first=is_first, qname=qname, flag=flag, rname=rname, mapq=mapq, cigar=cigar,
                                        rnext=rnext, tlen=tlen, seq=seq, qual=qual)

                elif line.startswith("@merged_read"):
                    if len(line[1:].split()) > 1:
                        merged_read_name, insert_str = line[1:].split()
                        inserts = [int(x) for x in insert_str.replace("inserts=", "").split(",")]
                    else:
                        merged_read_name = line[1:]
                        inserts = []

                    if merged_read_name != "merged_" + qname:
                        raise ValueError("Something wrong with the read names")

                    merged_read_seq = fh_in.next().rstrip()
                    sep = fh_in.next()
                    merged_read_qual = fh_in.next().rstrip()
                    merge_test.add_merged_read(merged_read_seq, inserts)
                    fh_out.writelines(merge_test.create_prelim_csv_lines())
                    merge_tests.append(merge_test)
                else:
                    raise ValueError("Shouldn't get here, invalid test fastq format. " + str(i) + ": " + line)

        return sam, merge_tests


    def test_get_sliced_merged_reads_from_prelim_csv(self):
        # Test insertions wrt ref.  No padding wrt ref.  No padding wrt slice.
        # Merge start = 4
        # Merge end = 13
        # "---GCCTAGACAT--------"  If there is no insertions.
        # 1-based positions
        MERGE_START = 4
        MERGE_END = 13
        SLICE_START_BEFORE_MERGE_START = MERGE_START - 1
        SLICE_START_AFTER_MERGE_START = MERGE_START +1
        SLICE_END_BEFORE_MERGE_END = MERGE_END - 1
        SLICE_END_AFTER_MERGE_END = MERGE_END + 1

        prelim_csv, merge_tests = self.__generate_merged_read_sam(TestSliceAligned.TEST_PAD_INSERT_FASTQ)
        for slice_start, slice_end in [(SLICE_START_BEFORE_MERGE_START, SLICE_END_BEFORE_MERGE_END),
                                       (SLICE_START_BEFORE_MERGE_START, SLICE_END_AFTER_MERGE_END),
                                       (SLICE_START_AFTER_MERGE_START, SLICE_END_BEFORE_MERGE_END),
                                       (SLICE_START_AFTER_MERGE_START, SLICE_END_AFTER_MERGE_END)]:
            LOGGER.info("\n\nRunning tests for insertions, slice=[{}, {}]".format(slice_start, slice_end))

            all_stats = sam.align_stats.AlignStats()

            sam_record = SamRecord()
            for i, (merged_name, merged_seq, read_stats) in enumerate(
                    pair.PairedRecord.create_msa_slice_from_sam(prelim_csv=prelim_csv, target_ref=TestSliceAligned.REFERENCE,
                                                           ref_len=TestSliceAligned.REF_LEN,
                                                           q_cutoff=TestSliceAligned.Q_CUTOFF,
                                                           mapping_cutoff=TestSliceAligned.M_CUTOFF,
                                                           pad_space_btn_mates="N",
                                                           is_pad_wrt_ref=False,
                                                           is_pad_wrt_slice=False,
                                                           is_insert=True, slice_start=slice_start,
                                                           slice_end=slice_end)):
                all_stats.merge_stats(read_stats)
                merge_test = merge_tests[i]
                sliced_merge_test_seq = merge_test.get_sliced_merged_read(slice_start, slice_end).lstrip("-").rstrip("-")
                self.assertEqual(merge_test.merged_read_name, merged_name,
                                 "Testcases out of order.  Expected testcase=" + merge_test.merged_read_name +
                                 " but got " + merged_name)

                LOGGER.debug( str(i) + ": " + merged_name + "=" + merged_seq)
                self.assertEquals(sliced_merge_test_seq, merged_seq,
                                  "Expected=" + sliced_merge_test_seq + " but " + merged_name + " actual=" + merged_seq)
            LOGGER.debug(all_stats.dump_stats())



    def test_get_merged_reads_from_prelim_csv(self):


        LOGGER.info( "Running tests for padded insertions merged sequence")
        prelim_csv, merge_tests = self.__generate_prelim_csv(TestSliceAligned.TEST_PAD_INSERT_FASTQ)
        #  Test merging reads with padding wrt ref and insertions wrt ref
        all_stats = sam.align_stats.AlignStats()
        for i, (merged_name, merged_seq, read_stats) in enumerate(
                slice.get_merged_reads_from_prelim_csv(prelim_csv=prelim_csv, target_ref=TestSliceAligned.REFERENCE,
                                                       ref_len=TestSliceAligned.REF_LEN,
                                                       q_cutoff=TestSliceAligned.Q_CUTOFF,
                                                       mapping_cutoff=TestSliceAligned.M_CUTOFF,
                                                       pad_space_btn_mates="N",
                                                       is_pad_wrt_ref=True,
                                                       is_insert=True,
                                                       is_pad_wrt_slice=False,
                                                       slice_start=None, slice_end=None)):
            merge_test = merge_tests[i]
            all_stats.merge_stats(read_stats)
            self.assertEqual(merge_test.merged_read_name, merged_name,
                             "Testcases out of order.  Expected testcase=" + merge_test.merged_read_name +
                             " but got " + merged_name)

            LOGGER.debug( str(i) + ": " + merged_name + "=" + merged_seq)
            self.assertEquals(merge_test.merged_read_seq, merged_seq,
                              "Expected=" + merge_test.merged_read_seq + " but " + merged_name + " actual=" + merged_seq)

        LOGGER.debug(all_stats.dump_stats())




if __name__ == '__main__':
    unittest.main()
