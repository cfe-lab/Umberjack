import unittest
import slice_miseq
import os
import sys
import csv
import array



SAM_FILENAME = "./simulations/data/indelible/cov1x/sample_genomes.100.1x.errfree.consensus.sam"
MSA_FASTA_FILENAME = "./simulations/out/consensus/mut100_1x_errfree_window400/sample_genomes.100.1x.consensus.consensus.msa.fasta"

MAPQ_CUTOFF = 20  # alignment quality cutoff
MAX_PROP_N = 0.1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/indelible/sample_genomes.100.consensus.fasta"
REF = "consensus"
REF_LEN = 9000
DNDS_DIR = "./simulations/out/consensus/mut100_1x_errfree_window400"
PVALUE = 0.05

MIN_WINDOW_BREADTH_COV_FRACTION = 0.875
MIN_WINDOW_DEPTH_COV = 50

WINDOWSIZE = 400
THREADS_PER_WINDOW = 2
WINDOW_PROCS = 2
START_NUCPOS = 3001
END_NUCPOS = 3400


SMOOTH_DIST = 10  # amino acid sites before or after the site to smooth over

ACTUAL_DNDS_FILENAME = DNDS_DIR + os.sep + 'actual_dnds_by_site.csv'
NUC_PER_CODON = 3

class TestSliceMiSeq(unittest.TestCase):




    def test_get_seq_dnds(self):
        with open(ACTUAL_DNDS_FILENAME, 'w') as dnds_fh:
            dnds_tsv_comments = ("ref=" + REF + "," +
                                 "ref_len=" + str(REF_LEN) + "," +
                                 "sam=" + SAM_FILENAME + "," +
                                 "mapping qual cutoff=" + str(MAPQ_CUTOFF) + "," +
                                 "read qual cutoff=" + str(READ_QUAL_CUTOFF) + "," +
                                 "max fraction N=" + str(MAX_PROP_N) + "," +
                                 "start nuc pos=" + str(START_NUCPOS) + "," +
                                 "end nuc pos=" + str(END_NUCPOS) + "," +
                                 "windowsize=" + str(WINDOWSIZE) + "," +
                                 "window depth thresh=" + str(MIN_WINDOW_DEPTH_COV) + "," +
                                 "window breadth fraction=" + str(MIN_WINDOW_BREADTH_COV_FRACTION) + "," +
                                 "pvalue=" + str(PVALUE) +
                                 "smooth_dist=" + str(SMOOTH_DIST))
            slice_miseq.tabulate_dnds(dnds_tsv_dir=DNDS_DIR, ref=REF, ref_nuc_len=REF_LEN, pvalue_thresh=PVALUE,
                                      output_csv_filename=ACTUAL_DNDS_FILENAME, smooth_dist=SMOOTH_DIST,
                                      comments=dnds_tsv_comments)

        # Check that the file exists
        self.assertTrue(os.path.exists(ACTUAL_DNDS_FILENAME) and os.path.getsize(ACTUAL_DNDS_FILENAME) > 0,
                        msg="Unable to write to dn/ds summary csv file")



    def tabulate_nuc_subst(self):
        slice_miseq.tabulate_nuc_subst(output_csv_filename="/home/thuy/gitrepo/MutationPatterns/out/140415_M01841_0059_000000000-A64EA/all_nuc_subst.csv",
                                       comments="",
                                       nucmodelfit_dir="/home/thuy/gitrepo/MutationPatterns/out/140415_M01841_0059_000000000-A64EA")
if __name__ == '__main__':
    unittest.main()
