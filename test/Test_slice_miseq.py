import unittest
import slice_miseq
import os
import sys
import csv
import array


# For now, use a *.remap.sam file (paired end reads aligned to a consensus sequence with indels removed).
SAM_FILENAME = "./simulations/data/sample_genomes.grinder-reads.rename.sam"
MSA_FASTA_FILENAME = "./simulations/data/out/sample_genomes.grinder-reads.rename.msa.fasta"
MAPQ_CUTOFF = 10  # alignment quality cutoff
MAX_PROP_N = 0.25  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 10   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/sample_genomes.consensus.fas"
REF = "consensus"
REF_LEN = 9000
DNDS_DIR = "./simulations/data/out" + os.sep + REF
PVALUE = 0.05
BREADTH_THRESH = 0.75

ACTUAL_DNDS_FILENAME = DNDS_DIR + os.sep + 'actual_dnds_by_site.tsv'
GAMMA_DNDS_LOOKUP_FILENAME = "./simulations/data/sample_genomes.rates"
GAMMA_DNDS_LOOKUP_COL_INTERVAL = "Interval"
GAMMA_DNDS_LOOKUP_COL_SITE = "Site"
GAMMA_DNDS_LOOKUP_COL_SCALE_FACTOR = "Scaling_factor"
GAMMA_DNDS_LOOKUP_COL_RATE_CLASS = "Rate_class"
GAMMA_DNDS_LOOKUP_COL_DNDS = "Omega"


NUC_PER_CODON = 3

class TestSliceMiSeq(unittest.TestCase):

    @staticmethod
    def expected_dnds(dnds_lookup_filename):
        """
        Look dn/ds from sample_genomes.rates file.
        :rtype array of float: array where each element is a site dn/ds value.
        :param indelible_rates_filename str: full file path to indelible rates output file.
        :param scaling.factor int: the amount the the branch lengths are scaled
        """

        with open(dnds_lookup_filename, 'r') as lookup_fh:
            site0based_to_dnds = array.array('f')
            for row in csv.DictReader(lookup_fh, delimiter=','):
                site = row[GAMMA_DNDS_LOOKUP_COL_SITE]
                dnds = row[GAMMA_DNDS_LOOKUP_COL_DNDS]
                site0based_to_dnds.append(float(dnds))

            return site0based_to_dnds


    # def test_get_best_window_size_from_depth(self):
    #     # TEST CASE:  test the algo for calculating best window size
    #     # Uneven coverage, depth file has missing base positions (i.e.  base positions with zero coverage)
    #     # Smaller window sizes are sufficient in some parts of the genome, but larger window sizes required for others.
    #     # Must return the largest sufficient window size.
    #     expected_windowsize = 9
    #     depth_filename = "./data/Fake.depth"
    #     cov_thresh = 50
    #     check_windowsize = slice_miseq.get_best_window_size_from_depthfile(depth_filename=depth_filename, cov_thresh=cov_thresh) # TODO: implement
    #     self.assertEqual(expected_windowsize, check_windowsize,
    #                      msg="Expected windowsize=" + str(expected_windowsize) + " but got " + str(check_windowsize))
    #
    #
    # def test_get_best_window_size(self):
    #     # TEST CASE:  test conversion form sam to bam to depth file, then obtain best window size
    #     expected_window_size = 6  # TODO:  get the actual expected num
    #     check_windowsize = slice_miseq.get_best_window_size_by_read_cov(sam_filename=SAM_FILENAME, ref_filename=REFERENCE_FASTA, cov_thresh=26)
    #     self.assertEqual(check_windowsize, expected_window_size, "Expected window size=" + str(expected_window_size) +
    #                     " but got " + str(check_windowsize))
    #
    #
    # def test_slice_msa_fasta(self):
    #     # TEST CASE:  test that the sliced file is written properly to file
    #     start_pos = 8598
    #     end_pos = 8897
    #     expected_total_seq = 639
    #     slice_fasta_filename = "./simulations/data/sample_genomes.grinder-reads.rename.msa." + str(start_pos) + "_" + str(end_pos) + ".fasta"
    #
    #     total_seq = slice_miseq.create_slice_msa_fasta(fasta_filename=MSA_FASTA_FILENAME, out_fasta_filename=slice_fasta_filename,
    #                                                    start_pos=start_pos, end_pos=end_pos, breadth_thresh=BREADTH_THRESH)
    #
    #     self.assertTrue(os.path.exists(slice_fasta_filename) and os.path.getsize(slice_fasta_filename) > 0,
    #                     msg="Unable to write to sliced multiple-sequence aligned fasta file")
    #     self.assertEqual(total_seq, expected_total_seq, msg="Expected total sequences written to file = " +
    #                                                         str(expected_total_seq) + " but got " + str(total_seq))
    #
    #     expected_header = ">40607"
    #     expected_sequence = "CCAAGGCCTAATTTCGACTAACTCCAGACACTTCTCTGGGTTANC-TCGCTGGNGGANGCG-TAGNGGANATNGGAGNAACGNCNTNTNGTAATCTGACNNNNGTNTACTTNATNGTTGGTNGATCGTGTGNGGAACAACAGTNCATTGGCGTCGTAAGACCAAGATCATTGAGNGGTAACGATATTACCGCCATACCGAATAAGACGATGACGGTAACAGGTGCCTCCCTGCCGTCCGTGGTCTTTCCTAACAGCCGCCTTTGCGCTGCCAC---------------------------"
    #     is_found_expected_header = False
    #
    #     with open(slice_fasta_filename, 'r') as slice_fasta_fh:
    #         for line in slice_fasta_fh:
    #             line = line.rstrip()
    #             if line.split()[0] == expected_header:
    #                 is_found_expected_header = True
    #                 check_seq = slice_fasta_fh.next().rstrip()
    #                 self.assertEqual(check_seq, expected_sequence,
    #                                  msg="Expected sequence " + expected_header + " to be " + expected_sequence + " but got " + check_seq)
    #
    #     self.assertTrue(is_found_expected_header, msg="Expected to find seq " + expected_header + " but didn't")
    #

    def test_get_seq_dnds(self):
        # TODO:  automate check output
        # TODO: check multiple ref contigs
        seq_dnds_info = slice_miseq.get_seq_dnds_info(dnds_tsv_dir=DNDS_DIR, pvalue_thresh=PVALUE, ref=REF, ref_codon_len=REF_LEN/NUC_PER_CODON)
        with open(ACTUAL_DNDS_FILENAME, 'w') as dnds_fh:
            dnds_fh.write("Ref\tSite\tdNdS\tWindows\tCodons\tNonSyn\tSyn\tSubst\n")
            expected_site_2_dnds = TestSliceMiSeq.expected_dnds(GAMMA_DNDS_LOOKUP_FILENAME)

            for site in range(1, seq_dnds_info.get_seq_len() + 1):  # 1based codon sites
                site_dnds = seq_dnds_info.get_site_ave_dnds(site_1based=site)
                window = seq_dnds_info.get_site_window_cov(site_1based=site)
                reads = seq_dnds_info.get_site_ave_read_cov(site_1based=site)
                nonsyn = seq_dnds_info.get_site_ave_nonsyn_subs(site_1based=site)
                syn = seq_dnds_info.get_site_ave_syn_subs(site_1based=site)
                subs = seq_dnds_info.get_site_ave_subs(site_1based=site)
                line = "\t".join((REF, str(site), str(site_dnds), str(window), str(reads),  str(nonsyn), str(syn), str(subs)))
                dnds_fh.write(line + "\n")


if __name__ == '__main__':
    unittest.main()
