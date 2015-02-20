import unittest
import os
import logging
import sys

from sam import sam_handler


SAM_FILENAME = "./simulations/data/indelible/cov1x/sample_genomes.100.1x.consensus.sam"
MSA_FASTA_FILENAME = "./simulations/data/out/consensus/mut100_1x_window400/sample_genomes.100.1x.consensus.consensus.msa.fasta"
MAPQ_CUTOFF = 20  # alignment quality cutoff
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/indelible/sample_genomes.100.consensus.fasta"
REF = "consensus"
REF_LEN = 9000


class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

    def test_get_msa_fasta_from_sam(self):
        # Convert the sam alignments to the given reference to a multiple sequence alignment fasta file
        sam_handler.create_msa_fasta_from_sam(sam_filename=SAM_FILENAME, ref=REF, ref_len=REF_LEN,
                                              out_fasta_filename=MSA_FASTA_FILENAME, mapping_cutoff=MAPQ_CUTOFF,
                                              read_qual_cutoff=READ_QUAL_CUTOFF, max_prop_N=MAX_PROP_N)
        self.assertTrue(os.path.exists(MSA_FASTA_FILENAME) and os.path.getsize(MSA_FASTA_FILENAME) > 0,
                        msg="Unable to write to MSA-aligned fasta file")

        expect_seq_header = ">otu8542_read72"
        expect_seq = "GTACTTGGTATCTGCTTACCAGAGAAGNCAGAGAGCATTCCNGACGTATATATCTGGATTTCATCTNAGCGAGCTNGTCGGGGGTGTCCACCTGACCGTAANAGGAGTGAACTCCCTATTCCNTCCNNCTACCGGATTCNCCTNCATCGATTTATGCCATTGGGGGGTACGATCCAGAACCCGGGTTATGGGTATGCTCTTGCGCGTACGATGTCAAATCAGCGCACATGCAGCTGTCGCACGGGCCCGGCGTCGTGCNGGACCGTGCTCTCGATCTAAAATCGTCATTGCTCACGCAAAGAGAAANACGGGTGTATGCATTTGTCTATTCTTGCATTCGATGCATTGGGGATTAGCTCCTTTNCTATGGCGAATACTTCGTGCCATCATAACCCCTTGCAG"
        with open(MSA_FASTA_FILENAME, 'r') as msa_fasta_fh:
            found_header = False
            for line in msa_fasta_fh:
                if line.rstrip() == expect_seq_header:
                    found_header = True
                    check_seq = msa_fasta_fh.next().rstrip()
                    is_found = check_seq.find(expect_seq)
                    self.assertGreaterEqual(is_found, 0,
                                            "Expected " + expect_seq_header + " contains sequence " + expect_seq + " but got " + check_seq)

            self.assertTrue(found_header, "Expected header " + expect_seq_header + " but not found")

    def test_get_ave_coverage_from_bam(self):

        # TEST CASE:  From sam, get a sorted indexed bam.  Get the average coverage in a specified window.
        bam_filename = sam_handler.sam_to_sort_bam(sam_filename=SAM_FILENAME, ref_filename=REFERENCE_FASTA)
        ave_cov_in_window = sam_handler.get_ave_coverage_from_bam(bam_filename=bam_filename, pos_start=1, pos_end=300, ref="consensus")

        expected_ave_cov = 5628
        check_ave_cov = int(round(ave_cov_in_window, 0))

        self.assertEqual(expected_ave_cov, check_ave_cov, msg="Expected ave cov=" + str(expected_ave_cov) + " but got " + str(check_ave_cov))

if __name__ == '__main__':
    unittest.main()

