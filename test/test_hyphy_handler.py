import unittest
import os
import shutil
import tempfile
import hyphy.hyphy_handler as hyphy


THREADS = 2

class TestHyphyHandler(unittest.TestCase):

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
        self.tmpfasta.write("ACGACG\n")
        self.tmpfasta.write(">read2\n")
        self.tmpfasta.write("ACTACG\n")
        self.tmpfasta.write(">read3\n")
        self.tmpfasta.write("AGGACG\n")
        self.tmpfasta.write(">read4\n")
        self.tmpfasta.write("TCGTCG\n")
        self.tmpfasta.write(">read5\n")
        self.tmpfasta.write("TCCTCC\n")
        self.tmpfasta.flush()  # flush python buffer
        os.fsync(self.tmpfasta.file.fileno())  # flush os buffer to disk
        self.tmpfasta.close()

        self.tmptree = tempfile.NamedTemporaryFile('w', dir=TEST_DIR, suffix=".tree", delete=False)
        self.tmptree.write("(read3:0.00058,read5:0.40084,(read2:0.00058,(read1:0.00058,read4:0.00057):0.14650):0.20946);")
        self.tmptree.flush()  # flush python buffer
        os.fsync(self.tmptree.file.fileno())  # flush os buffer to disk
        self.tmptree.close()



    def test_calc_dnds(self):
        """
        Test that we get output with the right number of sites from hyphy dN/dS calculations.
        Test that the output does not get overwritten.
        :return:
        """
        dnds = hyphy.calc_dnds(codon_fasta_filename=self.tmpfasta.name, tree_filename=self.tmptree.name, threads=2)
        self.assertTrue(dnds and os.path.exists(dnds), "HyPhy sitewise dN/dS file was not created")

        # Check number of lines in file = total sites + 1 (header)
        total_lines = 0
        expected_num_sites = 2
        with open(dnds, 'rU') as fh_in:
            for line in fh_in:
                total_lines += 1
        actual_num_sites = total_lines - 1
        self.assertEqual(actual_num_sites, expected_num_sites,
                        "Expected " + str(expected_num_sites) + " sites but got " + str(actual_num_sites))


        # Check that dN/dS tsv file doesn't get ovewritten.
        dnds_again = hyphy.calc_dnds(codon_fasta_filename=self.tmpfasta.name, tree_filename=self.tmptree.name, threads=2)

        modify_time = os.path.getmtime(dnds)
        modify_time_again = os.path.getmtime(dnds_again)
        self.assertEqual(modify_time, modify_time_again, "HyPhy dN/dS sitewise tsv file overwritten")



if __name__ == '__main__':
    unittest.main()
