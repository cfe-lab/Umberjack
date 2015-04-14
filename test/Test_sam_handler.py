import unittest
import sam.sam_handler
import os
import subprocess
import Utility
import shutil
import tempfile
import sam_test_case
import config.settings as settings

settings.setup_logging()

# Simulation Configs
SIM_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations"
SIM_DATA_FILENAME_PREFIX = "umberjack_unittest"
SIM_DATA_DIR = SIM_DIR + os.sep + "data" + os.sep + SIM_DATA_FILENAME_PREFIX


SIM_PIPELINE_PY = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "sim_pipeline.py"

# INDELible dN/dS values that INDELible is aiming to simulate
INDELIBLE_DNDS_FILENAME = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.rates.csv"

# Sliding Window configs
POPN_CONSENSUS_FASTA =  SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.consensus.fasta"
REF = "consensus"

# Sam file for ART reads alignments against consensus of INDELible population
SIM_SAM = SIM_DATA_DIR + os.sep + "mixed" + os.sep + "aln" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.reads.consensus.bwa.sort.query"

OUT_DIR = SIM_DIR + os.sep + "out" + SIM_DATA_FILENAME_PREFIX + os.sep + REF

MAPQ_CUTOFF = 20  # alignment quality cutoff
MAX_PROP_N = 0.5  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]

MIN_WINDOW_BREADTH_COV_FRACTION = 0.875

TEST_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "out"

TEST_MERGE_FASTQ = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.merge.fq"
TEST_PAIR_SELECTION_SAM = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.pairselection.sam"
TEST_PAIR_SELECTION_TARGET_REF = "targetref"
EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.pairselection.msa.fasta"

class TestSamHandler(unittest.TestCase):

    def setUp(self):

        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR)
            os.makedirs(TEST_DIR)
        self.merge_testcases = sam_test_case.SamTestCase.parse_fastq(TEST_MERGE_FASTQ)

        self.TMPSAM_REF1 = "ref1"
        self.TMPSAM_REF1_LEN = 4

        self.TMPSAM_REF2 = "ref2"
        self.TMPSAM_REF2_LEN = 400

        self.TMPSAM_REF_DNE = "!@#$"
        self.TMPSAM_REF_DNE_LEN = None

        # Unknown Sort Order Sam
        self.tmpsam_unsort = tempfile.NamedTemporaryFile('w', suffix=".nosort.sam", dir=TEST_DIR, delete=False)
        self.tmpsam_unsort.write("@HD\tVN:0.0\tSO:unknown\n")
        self.tmpsam_unsort.write("@PG\tID:fakeid\tPN:programname\tCL:fake command line\tDS:desc\tVN:version0\n")
        self.tmpsam_unsort.write("@CO\trandom comment1\n")
        self.tmpsam_unsort.write("@CO\trandom comment2\n")
        self.tmpsam_unsort.write("@SQ\tSN:{}\tLN:{}\tSP:fake species\n".format(self.TMPSAM_REF1, self.TMPSAM_REF1_LEN))
        self.tmpsam_unsort.write("@SQ\tSN:{}\tLN:{}\tSP:fake species\n".format(self.TMPSAM_REF2, self.TMPSAM_REF2_LEN))
        self.tmpsam_unsort.write("read1\t4\t*\t0\t49\t*\t*\t*\t0\tACGT\tHHHH")
        self.tmpsam_unsort.flush()  # flush python buffer
        os.fsync(self.tmpsam_unsort.file.fileno())  # flush os buffer to disk
        self.tmpsam_unsort.close()

        # Queryname Sort Order Sam
        self.tmpsam_query = tempfile.NamedTemporaryFile('w', suffix=".sort.query.sam", dir=TEST_DIR, delete=False)
        self.tmpsam_query.write("@HD\tVN:0.0\tSO:queryname\n")
        self.tmpsam_query.write("@PG\tID:fakeid\tPN:programname\tCL:fake command line\tDS:desc\tVN:version0\n")
        self.tmpsam_query.write("@CO\trandom comment1\n")
        self.tmpsam_query.write("@CO\trandom comment2\n")
        self.tmpsam_query.write("@SQ\tSN:ref1\tLN:4\tSP:fake species\n")
        self.tmpsam_query.write("@SQ\tSN:ref2\tLN:400\tSP:fake species\n")
        self.tmpsam_query.write("read1\t4\t*\t0\t49\t*\t*\t*\t0\tACGT\tHHHH")
        self.tmpsam_query.flush()  # flush python buffer
        os.fsync(self.tmpsam_query.file.fileno())  # flush os buffer to disk
        self.tmpsam_query.close()


        # No header sam
        self.tmpsam_noheader = tempfile.NamedTemporaryFile('w', suffix=".noheader.sam", dir=TEST_DIR, delete=False)
        self.tmpsam_noheader.write("read1\t4\t*\t0\t49\t*\t*\t*\t0\tACGT\tHHHH")
        self.tmpsam_noheader.flush()  # flush python buffer
        os.fsync(self.tmpsam_noheader.file.fileno())  # flush os buffer to disk
        self.tmpsam_noheader.close()


    #
    #
    # def setUp(self):
    #     """
    #     Generate simulated data for unit tests
    #     """
    #     # Specify sim_pipeline.py configurations into separate config file
    #     config_filename = SIM_DATA_DIR + os.sep + "umberjack_unittest.conf"
    #     subprocess.check_call(["python", SIM_PIPELINE_PY, config_filename])
    #
    #     if not os.path.exists(OUT_DIR):
    #         os.makedirs(OUT_DIR)

    def test_get_reflen(self):
        """
        Tests that header is parsed for ref len properly.
        """
        # Test returns reference length or None if not found
        actual_ref1_len = sam.sam_handler.get_reflen(sam_filename=self.tmpsam_unsort.name, ref=self.TMPSAM_REF1)
        self.assertEqual(self.TMPSAM_REF1_LEN, actual_ref1_len, "Expected {} but got {} for {}".format(self.TMPSAM_REF1_LEN, actual_ref1_len, self.TMPSAM_REF1))

        actual_ref2_len = sam.sam_handler.get_reflen(sam_filename=self.tmpsam_unsort.name, ref=self.TMPSAM_REF2)
        self.assertEqual(self.TMPSAM_REF2_LEN, actual_ref2_len, "Expected {} but got {} for {}".format(self.TMPSAM_REF2_LEN, actual_ref2_len, self.TMPSAM_REF2))

        actual_ref_dne_len = sam.sam_handler.get_reflen(sam_filename=self.tmpsam_unsort.name, ref=self.TMPSAM_REF_DNE)
        self.assertEqual(self.TMPSAM_REF_DNE_LEN, actual_ref_dne_len, "Expected {} but got {} for {}".format(self.TMPSAM_REF_DNE_LEN, actual_ref_dne_len, self.TMPSAM_REF_DNE))


        # Test that a sam with no header returns None for reference length
        actual_ref_dne_len = sam.sam_handler.get_reflen(sam_filename=self.tmpsam_noheader.name, ref=self.TMPSAM_REF_DNE)
        self.assertEqual(actual_ref_dne_len, self.TMPSAM_REF_DNE_LEN, "Expected {} but got {} for {}".format(self.TMPSAM_REF_DNE_LEN, actual_ref_dne_len, self.TMPSAM_REF_DNE))




    def test_is_query_sort(self):
        """
        Test that header parsed properly for sam order
        """
        actual_unknown_sort = sam.sam_handler.is_query_sort(sam_filename=self.tmpsam_unsort.name)
        self.assertEqual(False, actual_unknown_sort, "Expected False but got {} for is_query_sort of unknown sort order".format(actual_unknown_sort))

        actual_query_sort = sam.sam_handler.is_query_sort(sam_filename=self.tmpsam_query.name)
        self.assertEqual(True, actual_query_sort, "Expected True but got {} for is_query_sort of queryname sort order".format(actual_query_sort))

         # Test that a sam with no header returns False for query sorted
        actual_noheader = sam.sam_handler.is_query_sort(sam_filename=self.tmpsam_noheader.name)
        self.assertEqual(False, actual_noheader, "Expected False but got {} for is_query_sort of sam with no header".format(actual_noheader))


    @staticmethod
    def diff_fasta_line(fasta1, fasta2):
        """
        Checks if fasta contents are the same.  Order matters.  Header contents matter.
        Line endings don't matter.
        Assumes both fastas have each sequence on 1 line.
        :return:  newline separated concatenation of the line in fasta1 and fasta2 that doesn't match or None if all lines match
        :rtype: str
        """
        with open(fasta1, 'rU') as fh1, open(fasta2, 'rU') as fh2:
            i = -1
            for i, line1 in enumerate(fh1):
                try:
                    line2 = fh2.next()
                    if line1.rstrip() != line2.rstrip():
                        return "line " + str(i+1) + ":\n" + line1 + "\n" + line2
                except StopIteration:  # in case fh2 is already at eof
                    return "line " + str(i+1) + ":\n" + line1 + "\n<eof>"

            try:    # check if fh2 has more lines than fh1
                line2 = fh2.next()
                if line2:
                    i += 1
                    return "line " + str(i+1) + ":\n<eof>\n" + line2
            except StopIteration:
                pass  # in case fh2 is already at eof


        return None


    def test_create_msa_slice_from_sam_pair_selection(self):
        """
        Tests that the sam_handler.create_msa_slice_from_sam() is iterating through the records
        and selecting the correct records for pairing.

        - Test Missing Mates
        - Test Unmapped Mates
        - Test Mates Mapped to wrong ref
        - Test Pair Mapped to Wrong Ref
        - Test Low Map Qual Threshold:
        - Test Secondary, Chimeric Alignments:


        CIGAR:  Should be tested in sam_record
        - test H, S, X, =, M, P
        -
        """

        ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA = TEST_DIR + os.sep + os.path.basename(TEST_PAIR_SELECTION_SAM).replace(".sam", ".msa.fasta")

        # Test that the pairs are selected correctly.   We don't care about slices, breadth thresholds or N's or masking stop codons here.
        # But we do care about mapping quality and target references.
        actual_written = sam.sam_handler.create_msa_slice_from_sam(sam_filename=TEST_PAIR_SELECTION_SAM,
                                                                   ref=TEST_PAIR_SELECTION_TARGET_REF,
                                                                   out_fasta_filename=ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA,
                                                                   mapping_cutoff=MAPQ_CUTOFF,
                                                                   read_qual_cutoff=READ_QUAL_CUTOFF, max_prop_N=1.0,
                                                                   breadth_thresh=0, start_pos=None, end_pos=None,
                                                                   do_insert_wrt_ref=False, do_mask_stop_codon=False)


        self.assertTrue(os.path.exists(ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA) and os.path.getsize(ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA) > 0,
                        ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA + " doesn't exist or is empty")


        diff_line = TestSamHandler.diff_fasta_line(EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA, ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA)
        self.assertIsNone(diff_line,
                        "Expected full msa fasta " + EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA + " different than " +
                        ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA + ":\n"  + str(diff_line))

        expected_written = Utility.get_total_seq_from_fasta(EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA)
        self.assertEqual(expected_written, actual_written,
                         "Expect total written seq {} but got {} from {}".format(expected_written, actual_written, ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA))


    def test_create_msa_slice_from_sam_unsorted(self):
        """
        Tests that the sam_handler.create_msa_slice_from_sam() requires a queryname sorted sam
        """
        actual_tmpsam_unsort_msa_fasta = TEST_DIR + os.sep + os.path.basename(self.tmpsam_unsort.name).replace(".sam", ".msa.fasta")
        self.assertRaises(ValueError, sam.sam_handler.create_msa_slice_from_sam,
                          sam_filename=self.tmpsam_unsort.name,
                          ref=TEST_PAIR_SELECTION_TARGET_REF,
                          out_fasta_filename=actual_tmpsam_unsort_msa_fasta,
                          mapping_cutoff=MAPQ_CUTOFF,
                          read_qual_cutoff=READ_QUAL_CUTOFF,
                          max_prop_N=1.0,
                          breadth_thresh=0,
                          start_pos=None, end_pos=None,
                          do_insert_wrt_ref=False,
                          do_mask_stop_codon=False)


    def __write_sam_testcase(self, testcases, samfile):
        """
        Writes the sam records to file for the list of SamTestCase
        """
        # Assume each testcase has the same reference to length dict SamTestCase.ref2len
        with open(samfile, 'w') as fh_out:
            fh_out.write("@HD\tVN:0.0\tSO:queryname\n")
            for ref, ref_len in testcases[0].ref2len.iteritems():
                fh_out.write("@SQ\tSN:{}\tLN:{}\n".format(ref, ref_len))
            for testcase in testcases:
                lines = testcase.create_sam_lines()
                fh_out.write(lines)


    def test_create_msa_slice_from_sam_maxpropN(self):
        """
        Tests that the sam_handler.create_msa_slice_from_sam() filters out sequences that have too many N's or gaps
        """

        # We only care that the sequences are filtered by fraction of N's.
        # We don't care about breadth thresholds or slicing.
        ACTUAL_TEST_MERGE_FULL_MSA_FASTA = TEST_DIR + os.sep + os.path.basename(TEST_MERGE_FASTQ).replace(".fq", ".msa.fasta")
        TEST_MERGE_SAM = TEST_DIR + os.sep + os.path.basename(TEST_MERGE_FASTQ).replace(".fq", ".sam")
        self.__write_sam_testcase(self.merge_testcases, TEST_MERGE_SAM)

        self.assertTrue(os.path.exists(TEST_MERGE_SAM) and os.path.getsize(TEST_MERGE_SAM) > 0,
                        "Expected test case sam for merging records " + TEST_MERGE_SAM + " does not exist or is empty")

        actual_written = sam.sam_handler.create_msa_slice_from_sam(sam_filename=TEST_MERGE_SAM,
                                                                   ref=self.merge_testcases[0].target_ref,
                                                                   out_fasta_filename=ACTUAL_TEST_MERGE_FULL_MSA_FASTA,
                                                                   mapping_cutoff=MAPQ_CUTOFF,
                                                                   read_qual_cutoff=READ_QUAL_CUTOFF,
                                                                   max_prop_N=MAX_PROP_N, breadth_thresh=0,
                                                                   start_pos=0, end_pos=0, do_insert_wrt_ref=True,
                                                                   do_mask_stop_codon=True)

        self.assertTrue(os.path.exists(ACTUAL_TEST_MERGE_FULL_MSA_FASTA) and os.path.getsize(ACTUAL_TEST_MERGE_FULL_MSA_FASTA) > 0,
                        ACTUAL_TEST_MERGE_FULL_MSA_FASTA + " doesn't exist or is empty")

        actual_header2seq = Utility.get_seq_dict(ACTUAL_TEST_MERGE_FULL_MSA_FASTA)

        expected_written = 0
        for testcase in self.merge_testcases:
            expected_seq = testcase.get_sliced_merged_read(slice_start_pos_wrt_ref_1based=None, slice_end_pos_wrt_ref_1based=None,
                                                           do_pad_wrt_slice=True, do_insert_wrt_ref=True, do_mask_stop_codon=True)
            actual_seq = actual_header2seq.get(testcase.read_name, None)

            if expected_seq.count("N") / float(len(expected_seq)) > MAX_PROP_N:
                self.assertIsNone(actual_seq,
                                  "Expect read " + testcase.read_name + " should not be in " + ACTUAL_TEST_MERGE_FULL_MSA_FASTA)
            else:
                expected_written += 1
                self.assertEqual(expected_seq, actual_seq,
                                 "Expected {} but got {} for testcase {}".format(expected_seq, actual_seq, testcase.read_name))

        self.assertEqual(expected_written, actual_written,
                         "Expect total written seq {} but got {} from {}".format(expected_written, actual_written, ACTUAL_TEST_MERGE_FULL_MSA_FASTA))



if __name__ == '__main__':
    unittest.main()
