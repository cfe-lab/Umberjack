import unittest
import sam.sam_handler
import os
import subprocess
import Utility
import shutil
import tempfile

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


class TestSamHandler(unittest.TestCase):



    def setUp(self):
        """
        Generate simulated data for unit tests
        """
        # Specify sim_pipeline.py configurations into separate config file
        config_filename = SIM_DATA_DIR + os.sep + "slidingwindow_unittest.config"
        subprocess.check_call(["python", SIM_PIPELINE_PY, config_filename])

        if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)


    def test_create_msa_slice_from_sam_edge(self):
        """
        Edge cases for sam file.
        Empty sam file.
        Desired Reference not in sam.
        """
        pass


    def test_create_msa_slice_from_sam(self):
        mapping_cutoff = 20
        read_qual_cutoff = 20
        max_prop_N = 0.1
        breadth_thresh = 0.875
        start_pos = 2041
        end_pos = 2340
        is_insert = False
        ref_len = None
        ref = None

        START_POS = 30
        END_POS = 60
        OUT_SLICE_FASTA_FILENAME = OUT_DIR + os.sep + "{}{}_{}.fasta".format(SIM_DATA_FILENAME_PREFIX, START_POS, END_POS)

        written = sam.sam_handler.create_msa_slice_from_sam(sam_filename=SIM_SAM, ref=POPN_CONSENSUS_FASTA,
                                                            out_fasta_filename=OUT_SLICE_FASTA_FILENAME,
                                                            mapping_cutoff=MAPQ_CUTOFF,
                                                            read_qual_cutoff=READ_QUAL_CUTOFF,
                                                            max_prop_N=MAX_PROP_N, breadth_thresh=MIN_WINDOW_BREADTH_COV_FRACTION,
                                                            start_pos=START_POS, end_pos=END_POS,
                                                            is_insert=False)
        self.assertGreater(written, 0, "Expect lines written to " + OUT_SLICE_FASTA_FILENAME + " but none written")
        self.assertTrue(os.path.exists(OUT_SLICE_FASTA_FILENAME) and os.path.getsize(OUT_SLICE_FASTA_FILENAME) > 0,
                        OUT_SLICE_FASTA_FILENAME + " doesn't exist or is empty")


    def test_create_msa_slice_from_sam_merge(self):
        """
        Test that the paired-end sequence merge feature masks low quality bases, bases that conflict between mates.
        """

        samfile = tempfile.NamedTemporaryFile("w", delete=False)
        # @HD	VN:1.4	SO:queryname
        # @SQ	SN:consensus	LN:450
        from sam.sam_constants import SamFlag as SamFlag
        flag_1st_paired = SamFlag.IS_PAIRED | SamFlag.IS_MAPPED_IN_PROPER_PAIR | SamFlag.IS_FIRST | SamFlag.IS_MATE_REVERSE
        samfile.write("@HD\tVN:1.4\tSO:queryname]\n")
        samfile.write("@SQ\tSN:consensus\tLN:10\n")
        #ACGTACGTACGT
        #
        samfile.write("read1\t" + str(flag_1st_paired) + "\tconsensus\t1\t40\t" + \n")

if __name__ == '__main__':
    unittest.main()
