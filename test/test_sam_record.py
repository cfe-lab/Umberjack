import unittest
from sam.paired_records import PairedRecord
from sam.sam_record import SamRecord
import logging
import sys
import os
import sam_test_case
import config.settings as settings

LOGGER = logging.getLogger(__name__)


class TestSamRecord(unittest.TestCase):

    TEST_PAD_INSERT_FASTQ = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.pad.insert.fq"
    Q_CUTOFF = 20
    M_CUTOFF = 20


    def setUp(self):
        self.testcases = sam_test_case.SamTestCase.parse_fastq(TestSamRecord.TEST_PAD_INSERT_FASTQ)

    def test_slice_inserts_merged_reads(self):
        # Test insertions wrt ref, slice wrt ref

        for i, testcase in enumerate(self.testcases):
            LOGGER.debug("TestCase " + str(i) + ". Inserts. Actual: " + testcase.read_name)
            sam_rec1 = SamRecord(testcase.mate1.ref_len)
            sam_rec1.fill_record_vals(qname=testcase.mate1.qname, flag=testcase.mate1.flag, rname=testcase.mate1.rname,
                                      seq=testcase.mate1.seq, cigar=testcase.mate1.cigar, mapq=testcase.mate1.mapq,
                                      qual=testcase.mate1.qual, pos=testcase.mate1.pos,
                                      rnext=testcase.mate1.rnext, pnext=testcase.mate1.pnext)
            sam_rec2 = SamRecord(testcase.mate2.ref_len)
            sam_rec2.fill_record_vals(qname=testcase.mate2.qname, flag=testcase.mate2.flag, rname=testcase.mate2.rname,
                                      seq=testcase.mate2.seq, cigar=testcase.mate2.cigar, mapq=testcase.mate2.mapq,
                                      qual=testcase.mate2.qual, pos=testcase.mate2.pos,
                                      rnext=testcase.mate2.rnext, pnext=testcase.mate2.pnext)

            for do_pad_wrt_slice in [True, False]:
                if not testcase.is_same_ref() and testcase.mate1.rname == testcase.target_ref:
                    merged_seq, merged_qual, stats = sam_rec1.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=do_pad_wrt_slice,
                                                                           do_mask_low_qual=True, q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                           slice_start_wrt_ref_1based=testcase.slice_start,
                                                                           slice_end_wrt_ref_1based=testcase.slice_end,
                                                                           do_insert_wrt_ref=True)
                elif not testcase.is_same_ref() and testcase.mate2.rname == testcase.target_ref:
                    # Test no padding wrt to ref, no padding wrt slice, include inserts wrt ref
                    merged_seq, merged_qual, stats = sam_rec2.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=do_pad_wrt_slice,
                                                                           do_mask_low_qual=True, q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                           slice_start_wrt_ref_1based=testcase.slice_start,
                                                                           slice_end_wrt_ref_1based=testcase.slice_end,
                                                                           do_insert_wrt_ref=True)
                else:
                    paired_rec = PairedRecord(sam_rec1, sam_rec2)
                    merged_seq, stats = paired_rec.merge_sam_reads(q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                   pad_space_btn_mates="N",
                                                                   do_insert_wrt_ref=True, do_pad_wrt_ref=False, do_pad_wrt_slice=do_pad_wrt_slice,
                                                                   slice_start_wrt_ref_1based=testcase.slice_start, slice_end_wrt_ref_1based=testcase.slice_end)

                expected_merged_seq = testcase.get_sliced_merged_read(testcase.slice_start, testcase.slice_end, do_pad_wrt_slice=do_pad_wrt_slice)

                LOGGER.debug("TestCase " + str(i) + " - Inserts, Padding Wrt Slice=" + str(do_pad_wrt_slice) +
                             ". Actual: " + testcase.read_name + "=" + merged_seq + "\n" + stats.dump_stats())
                self.assertEquals(expected_merged_seq, merged_seq,
                                  "Expected=" + expected_merged_seq +
                                  " Actual= " + merged_seq +
                                  " MergeTest=" + testcase.test_desc)


    def test_slice_noinserts_merged_reads(self):
        # Test strip insertions wrt ref, slice wrt ref

        for i, testcase in enumerate(self.testcases):
            LOGGER.debug("TestCase " + str(i) + ". No Inserts. Actual: " + testcase.read_name)
            sam_rec1 = SamRecord(testcase.mate1.ref_len)
            sam_rec1.fill_record_vals(qname=testcase.mate1.qname, flag=testcase.mate1.flag, rname=testcase.mate1.rname,
                                      seq=testcase.mate1.seq, cigar=testcase.mate1.cigar, mapq=testcase.mate1.mapq,
                                      qual=testcase.mate1.qual, pos=testcase.mate1.pos,
                                      rnext=testcase.mate1.rnext, pnext=testcase.mate1.pnext)
            sam_rec2 = SamRecord(testcase.mate2.ref_len)
            sam_rec2.fill_record_vals(qname=testcase.mate2.qname, flag=testcase.mate2.flag, rname=testcase.mate2.rname,
                                      seq=testcase.mate2.seq, cigar=testcase.mate2.cigar, mapq=testcase.mate2.mapq,
                                      qual=testcase.mate2.qual, pos=testcase.mate2.pos,
                                      rnext=testcase.mate2.rnext, pnext=testcase.mate2.pnext)

            for do_pad_wrt_slice in [True, False]:
                if not testcase.is_same_ref() and testcase.mate1.rname == testcase.target_ref:
                    merged_seq, merged_qual, stats = sam_rec1.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=do_pad_wrt_slice,
                                                                           do_mask_low_qual=True, q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                           slice_start_wrt_ref_1based=testcase.slice_start,
                                                                           slice_end_wrt_ref_1based=testcase.slice_end,
                                                                           do_insert_wrt_ref=False)
                elif not testcase.is_same_ref() and testcase.mate2.rname == testcase.target_ref:
                    # Test no padding wrt to ref, no padding wrt slice, include inserts wrt ref
                    merged_seq, merged_qual, stats = sam_rec2.get_seq_qual(do_pad_wrt_ref=False, do_pad_wrt_slice=do_pad_wrt_slice,
                                                                           do_mask_low_qual=True, q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                           slice_start_wrt_ref_1based=testcase.slice_start,
                                                                           slice_end_wrt_ref_1based=testcase.slice_end,
                                                                           do_insert_wrt_ref=False)
                else:
                    paired_rec = PairedRecord(sam_rec1, sam_rec2)
                    merged_seq, stats = paired_rec.merge_sam_reads(q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                   pad_space_btn_mates="N",
                                                                   do_insert_wrt_ref=False, do_pad_wrt_ref=False, do_pad_wrt_slice=do_pad_wrt_slice,
                                                                   slice_start_wrt_ref_1based=testcase.slice_start, slice_end_wrt_ref_1based=testcase.slice_end)

                expected_merged_seq = testcase.get_sliced_merged_read(testcase.slice_start, testcase.slice_end,
                                                                      do_pad_wrt_slice=do_pad_wrt_slice, do_strip_inserts=True)

                LOGGER.debug("TestCase " + str(i) + " - No Inserts, Padding Wrt Slice=" + str(do_pad_wrt_slice) +
                             ". Actual: " + testcase.read_name + "=" + merged_seq + "\n" + stats.dump_stats())
                self.assertEquals(expected_merged_seq, merged_seq,
                                  "Expected=" + expected_merged_seq +
                                  " Actual= " + merged_seq +
                                  " MergeTest=" + testcase.test_desc)





if __name__ == '__main__':
    unittest.main()
