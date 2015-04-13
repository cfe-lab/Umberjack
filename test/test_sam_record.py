import unittest
from sam.paired_records import PairedRecord
from sam.single_record import SamRecord
import logging
import os
import sam_test_case
import config.settings as settings

settings.setup_logging()
LOGGER = logging.getLogger(__name__)

# TODO:  test stop codons with inserts
class TestSamRecord(unittest.TestCase):

    TEST_PAD_INSERT_FASTQ = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.merge.fq"
    Q_CUTOFF = 20
    M_CUTOFF = 20


    def setUp(self):
        self.testcases = sam_test_case.SamTestCase.parse_fastq(TestSamRecord.TEST_PAD_INSERT_FASTQ)

    def __slice_merged_reads(self, do_slice, do_pad_wrt_ref, do_pad_wrt_slice, do_mask_low_qual, do_insert_wrt_ref, do_mask_stop_codon):
        """
        Test merging sliced reads under various conditions.
        :param do_pad_wrt_ref:
        :param do_pad_wrt_slice:
        :param do_mask_low_qual:
        :param do_insert_wrt_ref:
        :param do_mask_stop_codon:
        :return:
        """



        for i, testcase in enumerate(self.testcases):
            if do_slice:
                slice_start = testcase.slice_start
                slice_end = testcase.slice_end
            else:
                slice_start = None
                slice_end = None

            conditions = (" do_pad_wrt_ref=" + str(do_pad_wrt_ref) +
                          " do_pad_wrt_slice=" + str(do_pad_wrt_slice) +
                          " do_mask_low_qual=" + str(do_mask_low_qual) +
                          " do_insert_wrt_ref=" + str(do_insert_wrt_ref) +
                          " do_mask_stop_codon=" + str(do_mask_stop_codon) +
                          " slice_start=" + str(slice_start) +
                          " slice_end=" + str(slice_end))

            LOGGER.debug("TestCase " + str(i+1) + conditions + " read=" + testcase.read_name)
            sam_rec1 = SamRecord(testcase.ref2len[testcase.target_ref])
            sam_rec1.fill_record_vals(qname=testcase.mate1.qname, flag=testcase.mate1.flag, rname=testcase.mate1.rname,
                                      seq=testcase.mate1.seq, cigar=testcase.mate1.cigar, mapq=testcase.mate1.mapq,
                                      qual=testcase.mate1.qual, pos=testcase.mate1.pos,
                                      rnext=testcase.mate1.rnext, pnext=testcase.mate1.pnext)
            sam_rec2 = SamRecord(testcase.ref2len[testcase.target_ref])
            sam_rec2.fill_record_vals(qname=testcase.mate2.qname, flag=testcase.mate2.flag, rname=testcase.mate2.rname,
                                      seq=testcase.mate2.seq, cigar=testcase.mate2.cigar, mapq=testcase.mate2.mapq,
                                      qual=testcase.mate2.qual, pos=testcase.mate2.pos,
                                      rnext=testcase.mate2.rnext, pnext=testcase.mate2.pnext)


            if not testcase.is_same_ref() and testcase.mate1.rname == testcase.target_ref:
                merged_seq, merged_qual, stats = sam_rec1.get_seq_qual(do_pad_wrt_ref=do_pad_wrt_ref,
                                                                       do_pad_wrt_slice=do_pad_wrt_slice,
                                                                       do_mask_low_qual=do_mask_low_qual,
                                                                       q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                       slice_start_wrt_ref_1based=slice_start,
                                                                       slice_end_wrt_ref_1based=slice_end,
                                                                       do_insert_wrt_ref=do_insert_wrt_ref,
                                                                       do_mask_stop_codon=do_mask_stop_codon)
            elif not testcase.is_same_ref() and testcase.mate2.rname == testcase.target_ref:
                # Test no padding wrt to ref, no padding wrt slice, include inserts wrt ref
                merged_seq, merged_qual, stats = sam_rec2.get_seq_qual(do_pad_wrt_ref=do_pad_wrt_ref,
                                                                       do_pad_wrt_slice=do_pad_wrt_slice,
                                                                       do_mask_low_qual=do_mask_low_qual,
                                                                       q_cutoff=TestSamRecord.Q_CUTOFF,
                                                                       slice_start_wrt_ref_1based=slice_start,
                                                                       slice_end_wrt_ref_1based=slice_end,
                                                                       do_insert_wrt_ref=do_insert_wrt_ref,
                                                                       do_mask_stop_codon=do_mask_stop_codon)
            else:
                paired_rec = PairedRecord(sam_rec1, sam_rec2)
                merged_seq, stats = paired_rec.get_seq_qual(q_cutoff=TestSamRecord.Q_CUTOFF,
                                                               pad_space_btn_mates="N",
                                                               do_insert_wrt_ref=do_insert_wrt_ref, do_pad_wrt_ref=do_pad_wrt_ref, do_pad_wrt_slice=do_pad_wrt_slice,
                                                               do_mask_stop_codon=do_mask_stop_codon,
                                                               slice_start_wrt_ref_1based=slice_start, slice_end_wrt_ref_1based=slice_end)

            expected_merged_seq = testcase.get_sliced_merged_read(slice_start_pos_wrt_ref_1based=slice_start,
                                                                  slice_end_pos_wrt_ref_1based=slice_end,
                                                                  do_pad_wrt_slice=do_pad_wrt_slice, do_insert_wrt_ref=do_insert_wrt_ref,
                                                                  do_mask_stop_codon=do_mask_stop_codon)

            LOGGER.debug("TestCase " + str(i+1) + conditions + " read=" + testcase.read_name +
                         ". Actual: " + testcase.read_name + "=" + merged_seq + "\n" + stats.dump_stats())
            self.assertEquals(expected_merged_seq, merged_seq,
                              "Expected=" + expected_merged_seq +
                              " Actual= " + merged_seq +
                              " MergeTest=" + testcase.test_desc + " read=" + testcase.read_name)


    def test_merged_reads(self):
        """
        Test the the reads are merged properly.  No slicing.
        TODO:  test when there are inserts before the slice that introduce or remove stop codons
        """
        self.__slice_merged_reads(do_slice=False,
                                  do_pad_wrt_ref=False,
                                  do_pad_wrt_slice=True,
                                  do_mask_low_qual=True,
                                  do_insert_wrt_ref=True,
                                  do_mask_stop_codon=True)

    def test_slice_merged_reads_inserts(self):
        """
        Test that insertions with respect to reference are correctly included, masked, or excluded in sliced merged reads.
        Don't mask stop codons.
        :return:
        """
        self.__slice_merged_reads(do_slice=True,
                                  do_pad_wrt_ref=False,
                                  do_pad_wrt_slice=True,
                                  do_mask_low_qual=True,
                                  do_insert_wrt_ref=True,
                                  do_mask_stop_codon=False)


    def test_slice_merged_reads_noinserts(self):
        """
        Test that all insertions with respect to reference are stripped correctly in sliced merged reads.
        Don't mask stop codons.
        """
        self.__slice_merged_reads(do_slice=True,
                                  do_pad_wrt_ref=False,
                                  do_pad_wrt_slice=True,
                                  do_mask_low_qual=True,
                                  do_insert_wrt_ref=False,
                                  do_mask_stop_codon=False)



    def test_slice_merged_reads_noinserts_mask_stop(self):
        """
        Test that stop codons are masked correctly when there are no inserts.
        """
        self.__slice_merged_reads(do_slice=True,
                                  do_pad_wrt_ref=False,
                                  do_pad_wrt_slice=True,
                                  do_mask_low_qual=True,
                                  do_insert_wrt_ref=True,
                                  do_mask_stop_codon=True)


    def test_slice_merged_reads_inserts_mask_stop(self):
        """
        Test that stop codons are masked correctly when there are inserts.
        TODO:  test when there are inserts before the slice that introduce or remove stop codons
        """
        self.__slice_merged_reads(do_slice=True,
                                  do_pad_wrt_ref=False,
                                  do_pad_wrt_slice=True,
                                  do_mask_low_qual=True,
                                  do_insert_wrt_ref=True,
                                  do_mask_stop_codon=True)

if __name__ == '__main__':
    unittest.main()
