import unittest
import sliding_window_tree
import os
import subprocess
import Utility
import shutil


# Simulation Configs
SIM_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations"
SIM_DATA_FILENAME_PREFIX = "slidingwindow_unittest"
SIM_DATA_DIR = SIM_DIR + os.sep + "data" + os.sep + SIM_DATA_FILENAME_PREFIX


SIM_PIPELINE_PY = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "sim_pipeline.py"

# INDELible dN/dS values that INDELible is aiming to simulate
INDELIBLE_DNDS_FILENAME = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.rates.csv"

# Full population dN/dS
EXPECTED_DNDS_FILENAME = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.dnds.tsv"

# Sliding Window configs
POPN_CONSENSUS_FASTA =  SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.consensus.fasta"
REF = "consensus"


MAPQ_CUTOFF = 20  # alignment quality cutoff
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]

MIN_WINDOW_BREADTH_COV_FRACTION = 0.875
MIN_WINDOW_DEPTH_COV = 10
WINDOW_SIZE = 300
WINDOW_SLIDE = 30
SMOOTH_DIST=10



THREADS_PER_WINDOW = 4
WINDOW_PROCS = 3



class TestSlidingWindowTree(unittest.TestCase):

    def setUp(self):
        """
        Generate simulated data for unit tests
        """
        # Specify sim_pipeline.py configurations into separate config file
        config_filename = SIM_DATA_DIR + os.sep + "slidingwindow_unittest.config"
        subprocess.check_call(["python", SIM_PIPELINE_PY, config_filename])



    def test_eval_windows_async(self):
        # ART generated reads aligned to population consensus
        SAM_FILENAME = SIM_DATA_DIR + os.sep + "mixed" + os.sep + "aln" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.reads.consensus.bwa.sort.query.sam"
        OUT_DIR =   SIM_DIR + os.sep + "out" + os.sep + SIM_DATA_FILENAME_PREFIX + os.sep + "Window" + str(WINDOW_SIZE)
        ACTUAL_DNDS_FILENAME = OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
        START_NUCPOS = 1
        END_NUCPOS = Utility.get_longest_seq_size_from_fasta(POPN_CONSENSUS_FASTA)
        # TODO:  automate check output of R scripts.  Right now, we need to manually view HTML generated from R.
        # i.e.  it's up to you to open up ./simulations/R/sliding_window_tree_unit_test.html and inspect the graphs/contents.
        sliding_window_tree.eval_windows_async(ref=REF, sam_filename=SAM_FILENAME,
                                               out_dir=OUT_DIR, map_qual_cutoff=MAPQ_CUTOFF,
                                               read_qual_cutoff=READ_QUAL_CUTOFF, max_prop_n=MAX_PROP_N,
                                               start_nucpos=START_NUCPOS, end_nucpos=END_NUCPOS,
                                               window_size=WINDOW_SIZE,
                                               window_depth_cutoff=MIN_WINDOW_DEPTH_COV,
                                               window_breadth_cutoff=MIN_WINDOW_BREADTH_COV_FRACTION,
                                               threads_per_window=THREADS_PER_WINDOW,
                                               concurrent_windows=WINDOW_PROCS,
                                               output_csv_filename=ACTUAL_DNDS_FILENAME,
                                               smooth_dist=SMOOTH_DIST, window_slide=WINDOW_SLIDE)

        rconfig_file = os.path.dirname(__file__) + os.sep +"simulations" + os.sep + "R" + os.sep + "sliding_window_tree_unit_test.config"
        with open(rconfig_file, 'w') as fh_out_config:
            fh_out_config.write("ACTUAL_DNDS_FILENAME=" + ACTUAL_DNDS_FILENAME + "\n")
            fh_out_config.write("EXPECTED_DNDS_FILENAME=" + EXPECTED_DNDS_FILENAME + "\n")
            fh_out_config.write("INDELIBLE_DNDS_FILENAME=" + INDELIBLE_DNDS_FILENAME + "\n")

        subprocess.check_call(["Rscript", "-e", "library(knitr); setwd('./simulations/R'); spin('sliding_window_tree_unit_test.R')"],
                              shell=False, env=os.environ)
        shutil.copy("./simulations/R/sliding_window_tree_unit_test.html",
                    OUT_DIR + os.sep + "sliding_window_tree_unit_test.html")


    def test_eval_windows_async_errfree(self):
        ERR_FREE_SAM_FILENAME = SIM_DATA_DIR + os.sep + "mixed" + os.sep + "aln" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.reads.errFree.consensus.bwa.sort.query.sam"
        ERR_FREE_OUT_DIR = SIM_DIR + os.sep + "out" + os.sep + SIM_DATA_FILENAME_PREFIX + os.sep + "Window" + str(WINDOW_SIZE) + ".errFree"
        ERR_FREE_ACTUAL_DNDS_CSV = ERR_FREE_OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
        START_NUCPOS = 1
        END_NUCPOS = Utility.get_longest_seq_size_from_fasta(POPN_CONSENSUS_FASTA)
        sliding_window_tree.eval_windows_async(ref=REF,
                                               sam_filename=ERR_FREE_SAM_FILENAME,
                                               out_dir=ERR_FREE_OUT_DIR,
                                               map_qual_cutoff=MAPQ_CUTOFF,
                                               read_qual_cutoff=READ_QUAL_CUTOFF,
                                               max_prop_n=MAX_PROP_N,
                                               start_nucpos=START_NUCPOS,
                                               end_nucpos=END_NUCPOS,
                                               window_size=WINDOW_SIZE,
                                               window_depth_cutoff=MIN_WINDOW_DEPTH_COV,
                                               window_breadth_cutoff=MIN_WINDOW_BREADTH_COV_FRACTION,
                                               threads_per_window=THREADS_PER_WINDOW,
                                               concurrent_windows=WINDOW_PROCS,
                                               output_csv_filename=ERR_FREE_ACTUAL_DNDS_CSV,
                                               mode=sliding_window_tree.MODE_DNDS,
                                               window_slide=WINDOW_SLIDE,
                                               smooth_dist=SMOOTH_DIST)
        # TODO:  verify concordance and correlation
        rconfig_file = os.path.dirname(__file__) + os.sep +"simulations" + os.sep + "R" + os.sep + "sliding_window_tree_unit_test.config"
        with open(rconfig_file, 'w') as fh_out_config:
            fh_out_config.write("ACTUAL_DNDS_FILENAME=" + ERR_FREE_ACTUAL_DNDS_CSV + "\n")
            fh_out_config.write("EXPECTED_DNDS_FILENAME=" + EXPECTED_DNDS_FILENAME + "\n")
            fh_out_config.write("INDELIBLE_DNDS_FILENAME=" + INDELIBLE_DNDS_FILENAME + "\n")

        subprocess.check_call(["Rscript", "-e", "library(knitr); setwd('./simulations/R'); spin('sliding_window_tree_unit_test.R')"],
                              shell=False, env=os.environ)
        shutil.copy("./simulations/R/sliding_window_tree_unit_test.html",
                    ERR_FREE_OUT_DIR + os.sep + "sliding_window_tree_unit_test.html")



if __name__ == '__main__':
    unittest.main()
