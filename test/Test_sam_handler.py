import unittest
import sam.sam_handler
import os
import subprocess
import Utility
import shutil
import tempfile
import sam_test_case

# Simulation Configs
SIM_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations"
SIM_DATA_FILENAME_PREFIX = "slidingwindow_unittest"
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
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]

MIN_WINDOW_BREADTH_COV_FRACTION = 0.875

TEST_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "out"

TEST_PAD_INSERT_FASTQ = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.pad.insert.fq"
TEST_PAIR_SELECTION_SAM = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.pairselection.sam"
TEST_PAIR_SELECTION_TARGET_REF = "targetref"
EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data" + os.sep + "test.pairselection.msa.fasta"

class TestSamHandler(unittest.TestCase):

    def setUp(self):

        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR)
            os.makedirs(TEST_DIR)
        self.merge_testcases = sam_test_case.SamTestCase.parse_fastq(TEST_PAD_INSERT_FASTQ)

        self.TMPSAM_REF1 = "ref1"
        self.TMPSAM_REF1_LEN = 4

        self.TMPSAM_REF2 = "ref2"
        self.TMPSAM_REF2_LEN = 400

        self.TMPSAM_REF_DNE = "!@#$"
        self.TMPSAM_REF_DNE_LEN = None

        # Unknown Sort Order Sam
        self.tmpsam_unsort = tempfile.NamedTemporaryFile('w', dir=TEST_DIR, delete=False)
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
        self.tmpsam_query = tempfile.NamedTemporaryFile('w', dir=TEST_DIR, delete=False)
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
        self.tmpsam_noheader = tempfile.NamedTemporaryFile('w', dir=TEST_DIR, delete=False)
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
    #     config_filename = SIM_DATA_DIR + os.sep + "slidingwindow_unittest.config"
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
    def is_same_fasta(fasta1, fasta2):
        """
        Checks if fasta contents are the same.  Order matters.  Header contents matter.
        Line endings don't matter.
        Assumes both fastas have each sequence on 1 line.
        :return:
        """
        with open(fasta1, 'rU') as fh1, open(fasta2, 'rU') as fh2:
            for line1 in fh1:
                line2 = fh2.next()
                if line1.rstrip() != line2.rstrip():
                    return False

        return True


    def test_create_msa_slice_from_sam_pair_selection(self):
        """
        Tests that the sam_handler.create_msa_slice_from_sam() is iterating through the records
        and selecting the correct records for pairing.

        Test Missing Mates:
        - paired mapped record followed by single mapped record.
            Expected: paired record is merged into 1 sequence into fasta.  single record has 1 sequence in fasta.
        - single mapped record followed by paired mapped record.
            Expected: paired record is merged into 1 sequence in fasta.  single record has 1 sequence in fasta.

        Test Unmapped Mates:
        - paired mapped record followed by paired record with only 1st mate mapped.
            Expected: paired record is merged into 1 sequence in fasta.  1st mapped mate in last pair has 1 sequence in fasta.
        - paired mapped record followed by paired record with only 2nd mate mapped.
            Expected: paired record is merged into 1 sequence in fasta.  2nd mapped mate in last pair has 1 sequence in fasta.
        - paired record with only 1st mate mapped followed by paired mapped record.
            Expected: paired record is merged into 1 sequence in fasta.  1st mapped mate in first pair has 1 sequence in fasta.
        - paired record with only 2nd mate mapped followed by paired mapped record.
            Expected: paired record is merged into 1 sequence in fasta.  2nd mapped mate in first pair has 1 sequence in fasta.


        Test Mates Mapped to wrong ref:
        - paired record with 1st mate mapped to target ref, 2nd mate mapped to wrong ref, followed by paired record mapped to right ref.
            Expected: paired record is merged into 1 sequence in fasta.  1st mapped mate in first pair has 1 sequence in fasta.
        - paired record with 1st mate mapped to wrong ref, 2nd mate mapped to right ref, followed by paired record mapped to right ref.
            Expected: paired record is merged into 1 sequence in fasta.  1st mapped mate in first pair has 1 sequence in fasta.

        Test Pair Mapped to Wrong Ref:
        - paired record mapped to wrong ref
            Expected: paired record not in fasta

        Test Low Map Qual Threshold:
        - paired record with 1st mate mapped to target ref with low qual, 2nd mate mapped to right ref with high qual, followed by paired record mapped to right ref with high qual.
            Expected: paired record is merged into 1 sequence in fasta.  1st mapped mate in first pair has 1 sequence in fasta.
        - paired record with 1st mate mapped to target ref with hi qual, 2nd mate mapped to right ref with low qual, followed by paired record mapped to right ref with high qual.
            Expected: paired record is merged into 1 sequence in fasta.  2nd mapped mate in first pair has 1 sequence in fasta.

        Test Multiple Alignments:
        - paired record - 1st mate secondary alignments on target ref & chimeric alignment on target ref, 2nd mate alignment on target ref, followed by
            paired record mapped to target ref
            Expected: 1st pair - 2nd mate sequence in fasta.  Last record single mate in fasta.
        - paired record - 2nd mate secondary alignments on target ref & chimeric alignment on target ref, 1st mate alignment on target ref, followed by
            paired record mapped to target ref
            Expected: 1st pair - 1st mate sequence in fasta.  Last record single mate in fasta.
        - paired record mapped to target ref followed by
            paired record - 1st mate secondary alignments on target ref,  with chimeric alignment on target ref, 2nd mate alignment on target ref
            Expected: 1st pair - 2nd mate sequence in fasta.  Last record single mate in fasta.
        - paired record mapped to target ref followed by
            paired record - 2nd mate secondary alignments on target ref with chimeric alignment on target ref, 1st mate alignment on target ref
            Expected: 1st pair - 1st mate sequence in fasta.  Last record single mate in fasta.


        CIGAR:  Should be tested in sam_record
        - test H, S, X, =, M, P
        -
        """

        START_POS = 30
        END_POS = 60
        ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA = TEST_DIR + os.sep + os.path.basename(TEST_PAIR_SELECTION_SAM).replace(".sam", ".msa.fasta")




        # Test that the pairs are selected correctly.   We don't care about slices, breadth thresholds or N's or masking stop codons here.
        # But we do care about mapping quality and target references.
        actual_written = sam.sam_handler.create_msa_slice_from_sam(sam_filename=TEST_PAIR_SELECTION_SAM,
                                                            ref=TEST_PAIR_SELECTION_TARGET_REF,
                                                            out_fasta_filename=ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA,
                                                            mapping_cutoff=MAPQ_CUTOFF,
                                                            read_qual_cutoff=READ_QUAL_CUTOFF,
                                                            max_prop_N=1.0,
                                                            breadth_thresh=0,
                                                            start_pos=None, end_pos=None,
                                                            is_insert=False,
                                                            is_mask_stop_codon=False)

        self.assertTrue(os.path.exists(ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA) and os.path.getsize(ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA) > 0,
                        ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA + " doesn't exist or is empty")

        expected_total_seq = Utility.get_total_seq_from_fasta(EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA)
        self.assertEqual(expected_total_seq, actual_written,
                           "Expected {} but got {} total sequences written to {}".format(expected_total_seq,
                                                                                         actual_written,
                                                                                         ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA))

        self.assertTrue(TestSamHandler.is_same_fasta(EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA, ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA),
                        "Expected full msa fasta " + EXPECTED_TEST_PAIR_SELECTION_FULL_MSA_FASTA + " different than " +
                        ACTUAL_TEST_PAIR_SELECTION_FULL_MSA_FASTA)


if __name__ == '__main__':
    unittest.main()
