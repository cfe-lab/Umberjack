import unittest
import unittest
import sliding_window_tree
import sys, os
import array
import csv
import subprocess
from hyphy import hyphy_handler as hyphy
from fasttree import fasttree_handler as fasttree
import Utility
import shutil

REFERENCE_FASTA =  os.path.dirname(__file__) + os.sep + "simulations/data/small/mixed/small.mixed.consensus.fasta"
REF = "consensus"
REF_LEN = Utility.get_seq2len(REFERENCE_FASTA)[REF]


SAM_FILENAME = os.path.dirname(__file__) + os.sep + "simulations/data/small/mixed/aln/small.mixed.reads.consensus.bowtie.sam"
MSA_FILENAME = SAM_FILENAME.replace(".sam", REF + ".msa.fasta")

MAPQ_CUTOFF = 20  # alignment quality cutoff
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]

MIN_WINDOW_BREADTH_COV_FRACTION = 0.875
MIN_WINDOW_DEPTH_COV = 10
WINDOW_SIZE = 300

# PVALUE = 0.05
# THREADS_PER_WINDOW = 4
# WINDOW_PROCS = 3
# START_NUCPOS = 3001
# END_NUCPOS = 3400
# SMOOTH_DIST=10
PVALUE = 0.05
THREADS_PER_WINDOW = 4
WINDOW_PROCS = 3
START_NUCPOS = 1
END_NUCPOS = REF_LEN
SMOOTH_DIST=10
WINDOW_SLIDE = 30


OUT_DIR =   os.path.dirname(__file__) + os.sep +"simulations/out/small/consensus/window" + str(WINDOW_SIZE)
ACTUAL_DNDS_FILENAME = OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
EXPECTED_DNDS_FILENAME = REFERENCE_FASTA.replace("consensus.fasta", "dnds.tsv")

HYPHY_BASEDIR = "/home/thuy/gitrepo/hyphy/res/TemplateBatchFiles"

class TestSlidingWindowTree(unittest.TestCase):

    def test_eval_windows_async(self):
        # TODO:  automate check output of R scripts.  Right now, we need to manually view HTML generated from R.
        # i.e.  it's up to you to open up ./simulations/R/sliding_window_tree_unit_test.html and inspect the graphs/contents.

        sliding_window_tree.eval_windows_async(ref=REF, sam_filename=SAM_FILENAME,
                                                               out_dir=OUT_DIR, map_qual_cutoff=MAPQ_CUTOFF,
                                                               read_qual_cutoff=READ_QUAL_CUTOFF, max_prop_n=MAX_PROP_N,
                                                               start_nucpos=START_NUCPOS, end_nucpos=END_NUCPOS,
                                                               window_size=WINDOW_SIZE,
                                                               window_depth_cutoff=MIN_WINDOW_DEPTH_COV,
                                                               window_breadth_cutoff=MIN_WINDOW_BREADTH_COV_FRACTION,
                                                               pvalue=PVALUE, threads_per_window=THREADS_PER_WINDOW,
                                                               concurrent_windows=WINDOW_PROCS,
                                                               output_csv_filename=ACTUAL_DNDS_FILENAME,
                                                               smooth_dist=SMOOTH_DIST, window_slide=WINDOW_SLIDE)

        rconfig_file = os.path.dirname(__file__) + os.sep +"simulations" + os.sep + "R" + os.sep + "sliding_window_unit_test.config"
        with open(rconfig_file, 'w') as fh_out_config:
            rconfig_file.write("ACTUAL_DNDS_FILENAME=" + ACTUAL_DNDS_FILENAME + "\n")
            rconfig_file.write("EXPECTED_DNDS_FILENAME=" + EXPECTED_DNDS_FILENAME + "\n")

        subprocess.check_call(["Rscript", "-e", "library(knitr); setwd('./simulations/R'); spin('sliding_window_tree_small_population_vhyphy.R')"],
                              shell=False, env=os.environ)
        shutil.copy("./simulations/R/sliding_window_tree_small_population_vhyphy.html",
                    OUT_DIR + os.sep + "sliding_window_tree_small_population_vhyphy.html")


    def test_eval_windows_async_errfree(self):
        ERR_FREE_ALN_CONSENSUS_SAM_FILENAME = SAM_FILENAME.replace(".reads.", ".reads.errFree.")
        ERR_FREE_OUT_DIR =   os.path.dirname(__file__) + os.sep +"simulations/out/small.errFree/consensus/window" + str(WINDOW_SIZE)
        ERR_FREE_OUTPUT_CSV = ERR_FREE_OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
        sliding_window_tree.eval_windows_async(ref=REF,
                                               sam_filename=ERR_FREE_ALN_CONSENSUS_SAM_FILENAME,
                                               out_dir=ERR_FREE_OUT_DIR,
                                               map_qual_cutoff=MAPQ_CUTOFF,
                                               read_qual_cutoff=READ_QUAL_CUTOFF,
                                               max_prop_n=MAX_PROP_N,
                                               start_nucpos=1,
                                               end_nucpos=END_NUCPOS,
                                               window_size=WINDOW_SIZE,
                                               window_depth_cutoff=MIN_WINDOW_DEPTH_COV,
                                               window_breadth_cutoff=MIN_WINDOW_BREADTH_COV_FRACTION,
                                               pvalue=PVALUE,
                                               threads_per_window=THREADS_PER_WINDOW,
                                               concurrent_windows=WINDOW_PROCS,
                                               output_csv_filename=ERR_FREE_OUTPUT_CSV,
                                               mode=sliding_window_tree.MODE_DNDS,
                                               window_slide=WINDOW_SLIDE,
                                               smooth_dist=SMOOTH_DIST)
        # TODO:  verify concordance and correlation
        rconfig_file = os.path.dirname(__file__) + os.sep +"simulations" + os.sep + "R" + os.sep + "sliding_window_unit_test.config"
        with open(rconfig_file, 'w') as fh_out_config:
            rconfig_file.write("ACTUAL_DNDS_FILENAME=" + ERR_FREE_OUTPUT_CSV + "\n")
            rconfig_file.write("EXPECTED_DNDS_FILENAME=" + EXPECTED_DNDS_FILENAME + "\n")

        subprocess.check_call(["Rscript", "-e", "library(knitr); setwd('./simulations/R'); spin('sliding_window_tree_small_population_vhyphy.R')"],
                              shell=False, env=os.environ)
        shutil.copy("./simulations/R/sliding_window_tree_small_population_vhyphy.html",
                    ERR_FREE_OUT_DIR + os.sep + "sliding_window_tree_small_population_vhyphy.html")


    # def test_full_popn_dnds(self):
    #     """
    #     Check that the dN/dS values in a window for a full population can be recovered using the window pipeline
    #      (fasta slice --> fasttree --> hyphy).
    #     We do not use reads from the population, we use the actual full population sequence.
    #     Repeat this for populations with various mutation rates.
    #     :return:
    #     """
    #     #MUTATION_SCALING_RATES = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    #     MUTATION_SCALING_RATES = [10.0]
    #     FULL_POPN_SEQ_FASTA_PATTERN = os.path.dirname(__file__) + "/simulations/out/checkFullPopulation/scaling_{}_TRUE.fasta"
    #     FULL_POPN_SEQ_SITE_DNDS_TABLE_PATTERN = os.path.dirname(__file__) + "/simulations/data/indelible/scaling_{}_RATES.txt"
    #
    #     FULL_POPN_START_NUCPOS = 1
    #     FULL_POPN_END_NUCPOS = 400
    #
    #     for scaling_rate in MUTATION_SCALING_RATES:
    #         full_popn_fasta = FULL_POPN_SEQ_FASTA_PATTERN.format(scaling_rate)
    #         full_popn_dnds = FULL_POPN_SEQ_SITE_DNDS_TABLE_PATTERN.format(scaling_rate)
    #
    #         sliding_window_tree.eval_window(msa_fasta_filename=full_popn_fasta, window_depth_cutoff=MIN_WINDOW_DEPTH_COV,
    #                                     window_breadth_cutoff=MIN_WINDOW_BREADTH_COV_FRACTION,
    #                                     start_window_nucpos=FULL_POPN_START_NUCPOS, end_window_nucpos=FULL_POPN_END_NUCPOS,
    #                                     pvalue=PVALUE, threads_per_window=THREADS_PER_WINDOW, mode="DNDS"
    #                                     #hyphy_exe=hyphy.HYPHY_EXE,
    #                                     #hyphy_basedir=HYPHY_BASEDIR,
    #                                     #fastree_exe=fasttree.FASTTREE_EXE
    #                                     )



if __name__ == '__main__':
    unittest.main()
