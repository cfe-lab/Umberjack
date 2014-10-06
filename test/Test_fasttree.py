import unittest
import os
import fasttree.fasttree_handler as fasttree
import constants

WINDOW_FASTA = "./simulations/data/out/consensus/mut100_1x_errfree_window400/sample_genomes.100.1x.errfree.consensus.consensus.msa.3001_3400.fasta"
WINDOW_TREE = "./simulations/data/out/consensus/mut100_1x_errfree_window400/sample_genomes.100.1x.errfree.consensus.consensus.msa.3001_3400.tree"
RATE_DIFF_THRESHOLD = 0.1

class MyTestCase(unittest.TestCase):
    def setUp(self):  # TODO:  create file if it doesn't exist
        self.assertTrue(os.path.exists(WINDOW_FASTA) and os.path.getsize(WINDOW_FASTA))

    def test_make_tree(self):
        if os.path.exists(WINDOW_TREE):
            os.remove(WINDOW_TREE)

        fasttree.make_tree(fasta_fname=WINDOW_FASTA, threads=constants.THREADS)

        self.assertTrue(os.path.exists(WINDOW_TREE))  # TODO:  check contents of tree

    def test_extract_gtr_rates(self):
        self.test_make_tree()

        expected_rates = (1.0857, 9.0008, 1.2009, 1.0288, 10.6586, 1.0000)
        rates = fasttree.extract_gtr_rates(WINDOW_FASTA)
        for i, rate in enumerate(rates):
            rate_diff = abs(rate - expected_rates[i])
            self.assertTrue(rate_diff < RATE_DIFF_THRESHOLD,
                            "Expected rate[{}]={} +/- {} but got {}".format(i, expected_rates[i], RATE_DIFF_THRESHOLD, rate))



if __name__ == '__main__':
    unittest.main()
