import unittest
import sliding_window_tree
import sys, os
import array
import csv
import subprocess


SAM_FILENAME = "./simulations/data/sample_genomes.grinder-reads.rename.sam"
MAPQ_CUTOFF = 10  # alignment quality cutoff
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 10   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/sample_genomes.consensus.fas"
REF = "consensus"
REF_LEN = 9000
OUT_DIR = "./simulations/data/out" + os.sep + REF
ACTUAL_DNDS_FILENAME = OUT_DIR + os.sep + 'ave_dnds_by_site.tsv'

GAMMA_DNDS_LOOKUP_FILENAME = "./simulations/data/sample_genomes.rates"
GAMMA_DNDS_LOOKUP_COL_INTERVAL = "Interval"
GAMMA_DNDS_LOOKUP_COL_SITE = "Site"
GAMMA_DNDS_LOOKUP_COL_SCALE_FACTOR = "Scaling_factor"
GAMMA_DNDS_LOOKUP_COL_RATE_CLASS = "Rate_class"
GAMMA_DNDS_LOOKUP_COL_DNDS = "Omega"


MIN_WINDOW_BREADTH_COV_FRACTION = 0.75
MIN_WINDOW_DEPTH_COV = 50

BWA_EXE = "/home/thuy/programs/bwa/bwa-0.7.8/bwa"  # TODO:  do not hardcode
BOWTIE_EXE = ""
BOWTIE_BUILD_EXE=""

PVALUE = 0.05
THREADS_PER_WINDOW = 3
WINDOW_PROCS = 5
START_NUCPOS = 1
END_NUCPOS = 9000


class TestSlidingWindowTree(unittest.TestCase):

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

    # def test_process_windows(self):
    #
    #     startpos = sys.argv[1]
    #     endpos = sys.argv[2]
    #     actual_dnds_filename = './simulations/data/actual_dnds.tsv'
    #     with open(actual_dnds_filename, 'w') as dnds_fh:
    #         dnds_fh.write("ref\tsite\tdnds")
    #         expected_site_2_dnds = TestSlidingWindowTree.expected_dnds(GAMMA_DNDS_LOOKUP_FILENAME)
    #
    #         ref2SeqDnDsInfo = sliding_window_tree.process_windows(sam_filename=SAM_FILENAME,
    #                                                             ref_fasta_filename=REFERENCE_FASTA,
    #                                                             out_dir=OUT_DIR,
    #                                                             mapping_cutoff=MAPQ_CUTOFF,
    #                                                             read_qual_cutoff=READ_QUAL_CUTOFF,
    #                                                             max_prop_N=MAX_PROP_N,
    #                                                             window_breadth_thresh=MIN_WINDOW_BREADTH_COV_FRACTION,
    #                                                             window_depth_thresh=MIN_WINDOW_DEPTH_COV,
    #                                                             pvalue=PVALUE,
    #                                                             threads=THREADS)
    #
    #         for ref in ref2SeqDnDsInfo:
    #             for site in range(ref2SeqDnDsInfo[ref].get_seq_len()):
    #                 site_dnds = ref2SeqDnDsInfo[ref].get_site_ave_dnds(site)


    def test_eval_windows_async(self):

        with open(ACTUAL_DNDS_FILENAME, 'w') as dnds_fh:
            dnds_fh.write("Ref\tSite\tdNdS")
            expected_site_2_dnds = TestSlidingWindowTree.expected_dnds(GAMMA_DNDS_LOOKUP_FILENAME)

            seq_dnds_info = sliding_window_tree.eval_windows_async(sam_filename=SAM_FILENAME,
                                                                   ref=REF,
                                                                   ref_len=REF_LEN,
                                                                   out_dir=OUT_DIR,
                                                                   mapping_cutoff=MAPQ_CUTOFF,
                                                                   read_qual_cutoff=READ_QUAL_CUTOFF,
                                                                   max_prop_N=MAX_PROP_N,
                                                                   windowsize=300,
                                                                   window_breadth_thresh=MIN_WINDOW_BREADTH_COV_FRACTION,
                                                                   window_depth_thresh=MIN_WINDOW_DEPTH_COV,
                                                                   start_nucpos=START_NUCPOS,
                                                                   end_nucpos=END_NUCPOS,
                                                                   pvalue=PVALUE,
                                                                   threads_per_window=THREADS_PER_WINDOW,
                                                                   concurrent_windows=WINDOW_PROCS)

            for site in range(seq_dnds_info.get_seq_len()):
                site_dnds = seq_dnds_info.get_site_ave_dnds(site)
                dnds_fh.write(REF + "\t" + str(site) + "\t" + str(site_dnds) + "\n")


if __name__ == '__main__':
    unittest.main()
