import unittest
import os
import fasttree.fasttree_handler as fasttree
import tempfile
import shutil


RATE_DIFF_THRESHOLD = 0.1
THREADS = 2


class TestFastTreeHandler(unittest.TestCase):

    def setUp(self):
        """
        Put temporary fasta file under ./out directory.
        Delete existing ./out directory and remake it so that we start fresh.
        """
        TEST_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "out"
        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR)
            os.makedirs(TEST_DIR)

        self.tmpfasta = tempfile.NamedTemporaryFile('w', dir=TEST_DIR, suffix=".fasta", delete=False)
        self.tmpfasta.write(">read1\n")
        self.tmpfasta.write("ACGTACGT\n")
        self.tmpfasta.write(">read2\n")
        self.tmpfasta.write("ACGTACGC\n")
        self.tmpfasta.write(">read3\n")
        self.tmpfasta.write("ACGGACGC\n")
        self.tmpfasta.write(">read4\n")
        self.tmpfasta.write("ANNTACGT\n")
        self.tmpfasta.write(">read5\n")
        self.tmpfasta.write("TTGGACGC\n")
        self.tmpfasta.flush()  # flush python buffer
        os.fsync(self.tmpfasta.file.fileno())  # flush os buffer to disk
        self.tmpfasta.close()


    def test_make_tree(self):
        """
        Tests that fasttree_handler.make_tree() actually creates a tree.
        :return:
        """
        # Test that tree exists
        treefilename = fasttree.make_tree(fasta_fname=self.tmpfasta.name, threads=THREADS)
        modify_time = os.path.getmtime(treefilename)
        self.assertTrue(treefilename and os.path.exists(treefilename), "Tree not created")

        # Test that it doesn't overwite the existing tree or create a new tree
        treefilename_again = fasttree.make_tree(fasta_fname=self.tmpfasta.name, threads=THREADS)
        modify_time_again = os.path.getmtime(treefilename)
        self.assertEqual(treefilename, treefilename_again, "Expect that the treefile should be named " + treefilename)
        self.assertEqual(modify_time, modify_time_again, "Expect that tree " + treefilename + " does not get overwritten")


    def test_make_tree_debug(self):
        """
        Tests that fasttree_handler.make_tree(debug=True) actually creates a tree.
        Tests that it actually produces debug logging.
        :return:
        """
        # Test that tree exists
        treefilename = fasttree.make_tree(fasta_fname=self.tmpfasta.name, threads=THREADS, debug=True)
        modify_time = os.path.getmtime(treefilename)
        self.assertTrue(treefilename and os.path.exists(treefilename), "Tree not created")
        stdouterr_file = treefilename.replace(".tree", ".fasttree.stdouterr.txt")
        self.assertTrue(os.path.exists(stdouterr_file),
                        "Expect debug fasttree_handler.make_tree() outputs fasttree stdout/stderr to " + stdouterr_file)

        # Test that it doesn't overwite the existing tree or create a new tree
        treefilename_again = fasttree.make_tree(fasta_fname=self.tmpfasta.name, threads=THREADS, debug=True)
        modify_time_again = os.path.getmtime(treefilename)
        self.assertEqual(treefilename, treefilename_again, "Expect that the treefile should be named " + treefilename)
        self.assertEqual(modify_time, modify_time_again, "Expect that tree " + treefilename + " does not get overwritten")


    def test_extract_gtr_rates(self):
        treefilename = fasttree.make_tree(fasta_fname=self.tmpfasta.name, threads=THREADS)
        fasttree_log = treefilename.replace(".tree", ".fasttree.log")
        expected_rates = (0.0319, 0.0319, 1.2101, 0.0319, 2.1046, 1.0000)
        rates = fasttree.extract_gtr_rates(fasttree_log)
        for i, rate in enumerate(rates):
            rate_diff = abs(rate - expected_rates[i])
            self.assertTrue(rate_diff < RATE_DIFF_THRESHOLD,
                            "Expected rate[{}]={} +/- {} but got {}".format(i, expected_rates[i], RATE_DIFF_THRESHOLD, rate))



    def test_make_tree_repro(self):
        SINGLE_THREAD_FASTTREE_EXE = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "bin" + os.sep + "fasttree" + os.sep + "FastTree"
        treefilename = fasttree.make_tree(fasta_fname=self.tmpfasta.name, threads=THREADS)
        repro_treefilename = fasttree.make_tree_repro(fasta_fname=self.tmpfasta.name, intree_fname=treefilename, fastree_exe=SINGLE_THREAD_FASTTREE_EXE)

        # Test that repro tree exists
        self.assertTrue(repro_treefilename and os.path.exists(repro_treefilename), "Repro Tree not created")

        # Test that the input and output trees are the same
        with open(treefilename, 'rU') as fh_tree, open(repro_treefilename, 'rU') as fh_repro_tree:
            tree = fh_tree.readline()
            repro_tree = fh_repro_tree.readline()
            self.assertEqual(tree, repro_tree)


if __name__ == '__main__':
    unittest.main()
