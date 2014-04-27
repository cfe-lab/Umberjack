import unittest
import sam_handler
import os
import logging
import sys



SAM_FILENAME = "./simulations/data/sample_genomes.grinder-reads.rename.sam"
MSA_FASTA_FILENAME = "./simulations/data/out/consensus/sample_genomes.grinder-reads.rename.msa.fasta"
EXPECTED_MSA_FASTA_FILENAME = './simulations/data/sample_genomes.grinder-reads.rename.msa.expected.fasta'
MAPQ_CUTOFF = 0  # alignment quality cutoff
MAX_PROP_N = 0.2  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/sample_genomes.consensus.fas"
REF = "consensus"
REF_LEN = 9000

class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

    def test_get_msa_fasta_from_sam(self):
        # Convert the sam alignments to the given reference to a multiple sequence alignment fasta file
        sam_handler.create_msa_fasta_from_sam(sam_filename=SAM_FILENAME, ref=REF, ref_len=REF_LEN,
                                              out_fasta_filename=MSA_FASTA_FILENAME, mapping_cutoff=MAPQ_CUTOFF,
                                              read_qual_cutoff=READ_QUAL_CUTOFF, max_prop_N=MAX_PROP_N)
        self.assertTrue(os.path.exists(MSA_FASTA_FILENAME) and os.path.getsize(MSA_FASTA_FILENAME) > 0,
                        msg="Unable to write to MSA-aligned fasta file")


        # TEST CASE: mate1 and mate2 of read both hit the reference for the same length.  No padding required.
        #   Reverse complement of mate1 hits reference at base1.  Mate2 hits reference at base1.
        #

        expect_seq_header = ""
        expect_seq = ""
        with open(MSA_FASTA_FILENAME, 'r') as msa_fasta_fh:
            found_header = False
            for line in msa_fasta_fh:
                if line.rstrip() == expect_seq_header:
                    found_header = True
                    check_seq = msa_fasta_fh.next().rstrip()
                    self.assertEqual(check_seq, expect_seq,
                                     "Expected " + expect_seq_header + " has sequence " + expect_seq + " but got " + check_seq)


            self.assertTrue(found_header, "Expected header " + expect_seq_header + " but not found")

    def test_get_ave_coverage_from_bam(self):

        # TEST CASE:  From sam, get a sorted indexed bam.  Get the average coverage in a specified window.
        bam_filename = sam_handler.sam_to_sort_bam(sam_filename=SAM_FILENAME, ref_filename=REFERENCE_FASTA)
        ave_cov_in_window = sam_handler.get_ave_coverage_from_bam(bam_filename=bam_filename, pos_start=1, pos_end=300, ref="TestSample-RT")

        expected_ave_cov = 21
        check_ave_cov = int(round(ave_cov_in_window, 0))

        self.assertEqual(expected_ave_cov, check_ave_cov, msg="Expected ave cov=" + str(expected_ave_cov) + " but got " + str(check_ave_cov))
        # Only checked that the coverage as reported in the depth file matched what was output by the expected average coverage.
        # awk '{if ($2 <= 300) {lines++; total += $3}} END {ave=total/lines; print total " "  lines " " ave}  './data/TestSample-RT_S17.HIV1B-vif.remap.bam.sort.depth
        # 6237 300 20.79

if __name__ == '__main__':
    unittest.main()

