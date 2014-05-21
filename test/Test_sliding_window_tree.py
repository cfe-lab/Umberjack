import unittest
import sliding_window_tree
import sys, os
import array
import csv
import subprocess


SAM_FILENAME = "./simulations/data/indelible/cov1x/sample_genomes.100.1x.errfree.consensus.sam"
MAPQ_CUTOFF = 20  # alignment quality cutoff
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/indelible/sample_genomes.100.consensus.fasta"
REF = "consensus"
REF_LEN = 9000
OUT_DIR =  "./simulations/data/out/consensus/mut100_1x_errfree_window400"
ACTUAL_DNDS_FILENAME = OUT_DIR + os.sep + 'actual_dnds_by_site.tsv'

MIN_WINDOW_BREADTH_COV_FRACTION = 0.875
MIN_WINDOW_DEPTH_COV = 50
WINDOW_SIZE = 400

PVALUE = 0.05
THREADS_PER_WINDOW = 4
WINDOW_PROCS = 4
START_NUCPOS = 3001
END_NUCPOS = 3004



class TestSlidingWindowTree(unittest.TestCase):

    def test_eval_windows_async(self):
        # TODO:  automate check output of R scripts.  Right now, we need to manually view HTML generated from R.
        # i.e.  it's up to you to open up ./simulations/R/sliding_window_tree_unit_test.html and inspect the graphs/contents.

        seq_dnds_info = sliding_window_tree.eval_windows_async(sam_filename=SAM_FILENAME,
                                                               ref=REF,
                                                               ref_len=REF_LEN,
                                                               out_dir=OUT_DIR,
                                                               mapping_cutoff=MAPQ_CUTOFF,
                                                               read_qual_cutoff=READ_QUAL_CUTOFF,
                                                               max_prop_N=MAX_PROP_N,
                                                               windowsize=WINDOW_SIZE,
                                                               window_breadth_thresh=MIN_WINDOW_BREADTH_COV_FRACTION,
                                                               window_depth_thresh=MIN_WINDOW_DEPTH_COV,
                                                               start_nucpos=START_NUCPOS,
                                                               end_nucpos=END_NUCPOS,
                                                               pvalue=PVALUE,
                                                               threads_per_window=THREADS_PER_WINDOW,
                                                               concurrent_windows=WINDOW_PROCS,
                                                               output_dnds_tsv_filename=ACTUAL_DNDS_FILENAME)

        subprocess.check_call(["Rscript", "-e", "library(knitr); setwd('./simulations/R'); getwd(); spin('sliding_window_tree_unit_test.R')"],
                              shell=False, env=os.environ)






if __name__ == '__main__':
    unittest.main()
