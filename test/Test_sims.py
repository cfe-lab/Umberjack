# Tests that the simulated data are actually simulating what they're supposed to

# calculates the subsitution rate at the nucleotide level and amino acid level
# checks that the expected substitution rate is the actual substitution rate
import unittest
from Bio import Phylo
from Bio import Seq as Seq
import simulations.indelible.indelible_handler as indelibler
import simulations.phylo_dnds as phylo_dnds
from cStringIO import StringIO
import csv
import os
import subprocess


class TestSims(unittest.TestCase):

    @staticmethod
    def get_tree_len(treefilename):
        """
        Gets the tree length, i.e.  the sum of every branch in the tree.  Excludes root branch length.
        """
        tree = Phylo.read(treefilename, "newick")
        root_branch_length = tree.clade.branch_length
        tree_branch_length  = tree.clade.total_branch_length()
        return tree_branch_length - root_branch_length

    @staticmethod
    def print_tree_stats(treefilename):
        tree = Phylo.read(treefilename, "newick")

        clade_depths = tree.depths()
        root_branch_length = tree.clade.branch_length
        tree_branch_length  = tree.clade.total_branch_length()
        longest_depth = 0.0
        longest_tip = ""
        total_length= 0.0

        for clade, depth in clade_depths.iteritems():
            total_length += clade.branch_length
            if clade.is_terminal():
                print "leaf " + clade.name + " depth=" + str(depth)  + " distfromroot=" + str(tree.distance(clade)) + " branchlen=" + str(clade.branch_length)
                if longest_depth < depth:
                    longest_depth = depth
                    longest_tip = clade.name
            else:
                print "inner node " + str(clade.name) + " depth=" + str(depth)  + " distfromroot=" + str(tree.distance(clade)) + " branchlen=" + str(clade.branch_length)


        # NB:  the root node has depth=1.0.  Subtract that from the numbers so we're not confused.
        longest_depth -= root_branch_length
        total_length -= root_branch_length
        tree_branch_length -= root_branch_length
        print "longest depth = " + str(longest_depth) + " from " + longest_tip
        print "total length=" + str(total_length)  + " root_branch_length=" + str(root_branch_length) + " tree_branch_length=" + str(tree_branch_length)


    def test_indelible_tree_len(self):
        # TODO:  explicitly create indelible 50 scaling tree
        expected_tree_len = 50
        indelible_tree_txt = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small/trees.txt"
        indelible_tree_str = indelibler.get_tree_string(indelible_tree_txt)

        indelible_tree_handle = StringIO(indelible_tree_str)
        actual_tree_len = TestSims.get_tree_len(indelible_tree_handle)
        indelible_tree_handle.close()

        self.assertAlmostEqual(expected_tree_len, actual_tree_len, 0,
                               "Expected tree len=" + str(expected_tree_len) + " but got " + str(actual_tree_len))


    def test_calc_total_poss_subst(self):
        # ACT-T
        # A mutations: CCT-P, GCT-A, TCT-C.  0 syn, 3 nonsyn
        # C mutations: AAT-A, AGT-S, ATT-I.  0 syn, 3 nonsyn
        # T mutations: ACC-T, ACG-T, ACA-T.  3 syn, 0 nonsyn
        codon = Seq.Seq("ACT")
        expected_poss_syn = 3
        expected_poss_nonsyn = 6
        actual_poss_syn, actual_poss_nonsyn = phylo_dnds.calc_total_poss_subst(codon)
        self.assertEqual(actual_poss_syn, expected_poss_syn,
                         "Expected possible synonymous substitutions=" + str(expected_poss_syn) + " but got " + str(actual_poss_syn))
        self.assertEqual(actual_poss_nonsyn, expected_poss_nonsyn,
                         "Expected possible synonymous substitutions=" + str(expected_poss_nonsyn) + " but got " + str(actual_poss_nonsyn))


    def test_calc_subst(self):
        # ACT-T
        # A mutations: CCT-P, GCT-A, TCT-C.  0 syn, 3 nonsyn
        # C mutations: AAT-A, AGT-S, ATT-I.  0 syn, 3 nonsyn
        # T mutations: ACC-T, ACG-T, ACA-T.  3 syn, 0 nonsyn
        codon1 = Seq.Seq("ACT")
        codon2 = Seq.Seq("CCG")
        expected_syn = 1
        expected_nonsyn = 1
        total_syn, total_nonsyn = phylo_dnds.calc_total_subst(codon1, codon2)
        self.assertEqual(total_syn, expected_syn,
                         "Expected possible synonymous substitutions=" + str(expected_syn) + " but got " + str(total_syn))
        self.assertEqual(total_nonsyn, expected_nonsyn,
                         "Expected possible synonymous substitutions=" + str(expected_nonsyn) + " but got " + str(total_nonsyn))


    def test_calc_popn_dnds_simple(self):
        # GTA-V
        # G: ATA-I, TTA-L, CTA-L.  0 syn, 3 nonsyn
        # T: GAA-G, GGA-G, GCA-A.  0 syn, 3 nonsyn
        # A: GTG-V, GTC-V, GTT-V.  3 syn, 0 nonsyn
        ancestor_fasta_str = (">ROOT\nGTA\n" +
                              ">ANC1\nGTA\n")
        leaf_fasta_str =     (">LEF1\nGTG\n" +  # 1 syn
                              ">LEF2\nATA\n")  # 1 nonsyn
        tree_str = "((LEF1,LEF2)ANC1)ROOT;"

        expected_nonsyn = 1.0
        expected_syn = 1.0
        expected_poss_syn = 3.0
        expected_poss_nonsyn = 6.0
        expected_dnds = (expected_nonsyn/expected_poss_nonsyn)/(expected_syn/expected_poss_syn)

        leaf_handle = StringIO(leaf_fasta_str)
        ancestor_handle = StringIO(ancestor_fasta_str)
        tree_handle = StringIO(tree_str)
        actual_sites_dnds = phylo_dnds.calc_popn_dnds(leaf_fasta=leaf_handle, ancestor_fasta=ancestor_handle, treefile=tree_handle)
        leaf_handle.close()
        ancestor_handle.close()
        tree_handle.close()

        self.assertAlmostEqual(expected_dnds, actual_sites_dnds[0], 2,
                               "Site: 0" " expected=" + str(expected_dnds) + " actual=" + str(actual_sites_dnds[0]))


    def test_calc_popn_dnds_simtest(self):
        """
        Tests that the simulated datasets are generated properly.
        Small enough that dn/ds can be calculated manually.
        :return:
        """
        leaf_fasta = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/simtest/50.0/simtest.50.0_TRUE.fasta")
        ancestor_fasta = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/simtest/50.0/simtest.50.0_ANCESTRAL.fasta")
        expected_sites_dnds_csv = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/simtest/50.0/simtest.50.0_RATES.csv")
        indelible_tree_txt = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/simtest/50.0/trees.txt")

        # Generate the simulated files
        simtest_unit_test_config_file =  os.path.abspath(os.path.dirname(__file__) + "/simulations/data/simtest/simtest.config")
        sim_pipeline_exe =  os.path.abspath(os.path.dirname(__file__) + "/simulations/sim_pipeline.py")
        subprocess.check_call(["python", sim_pipeline_exe, simtest_unit_test_config_file], env=os.environ)
        self.assertTrue(os.path.exists(leaf_fasta))
        self.assertTrue(os.path.exists(ancestor_fasta))
        self.assertTrue(os.path.exists(expected_sites_dnds_csv))
        self.assertTrue(os.path.exists(indelible_tree_txt))


        tree_handle = indelibler.get_tree_stringio(indelible_tree_txt)
        actual_sites_dnds = phylo_dnds.calc_popn_dnds(leaf_fasta=leaf_fasta, ancestor_fasta=ancestor_fasta, treefile=tree_handle)
        tree_handle.close()

        with open(expected_sites_dnds_csv, 'rU') as fh_expected:
            expected_reader = csv.DictReader(fh_expected)  # Site	Class	Partition	Inserted?	Omega
            for expected_row in  expected_reader:
                site_1based = int(expected_row["Site"])
                expected_dnds = float(expected_row["Omega"])
                actual_dnds = actual_sites_dnds[site_1based-1]
                if round(expected_dnds, 2) != round(actual_dnds, 2):
                    print "Site: " + str(site_1based) + " ! expected=" + str(expected_dnds) + " actual=" + str(actual_dnds)
                else:
                    print "Site: " + str(site_1based) + " K expected=" + str(expected_dnds) + " actual=" + str(actual_dnds)

        with open(expected_sites_dnds_csv, 'rU') as fh_expected:
            expected_reader = csv.DictReader(fh_expected)  # Site	Class	Partition	Inserted?	Omega
            for expected_row in  expected_reader:
                site_1based = int(expected_row["Site"])
                expected_dnds = float(expected_row["Omega"])
                actual_dnds = actual_sites_dnds[site_1based-1]
                self.assertAlmostEqual(expected_dnds, actual_dnds, 2,
                                       "Site: " + str(site_1based) + " expected=" + str(expected_dnds) + " actual=" + str(actual_dnds))





    def test_calc_popn_dnds(self):
        """
        Generates the simulated unit tests for SlidingWindow.
        Checks that the dN/dS given by Indelible are actually generated by Indelible.
        Large enough that dN/dS can't be calculated manually and
            large enough that the SlidingWindow will have sufficient sequences to process each window.
        :return:
        """
        leaf_fasta = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/small/50.0/small.50.0_TRUE.fasta")
        ancestor_fasta = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/small/50.0/small.50.0_ANCESTRAL.fasta")
        expected_sites_dnds_csv = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/small/50.0/small.50.0_RATES.csv")
        indelible_tree_txt = os.path.abspath(os.path.dirname(__file__) + "/simulations/data/small/trees.txt")

        # Generate the simulated files
        small_unit_test_config_file =  os.path.abspath(os.path.dirname(__file__) + "/simulations/data/small/small.config")
        sim_pipeline_exe =  os.path.abspath(os.path.dirname(__file__) + "/simulations/sim_pipeline.py")
        subprocess.check_call(["python", sim_pipeline_exe, small_unit_test_config_file], env=os.environ)
        self.assertTrue(os.path.exists(leaf_fasta))
        self.assertTrue(os.path.exists(ancestor_fasta))
        self.assertTrue(os.path.exists(expected_sites_dnds_csv))
        self.assertTrue(os.path.exists(indelible_tree_txt))


        indelible_tree_handle = indelibler.get_tree_stringio(indelible_tree_txt)
        actual_sites_dnds = phylo_dnds.calc_popn_dnds(leaf_fasta=leaf_fasta, ancestor_fasta=ancestor_fasta, treefile=indelible_tree_handle)
        indelible_tree_handle.close()

        with open(expected_sites_dnds_csv, 'rU') as fh_expected:
            expected_reader = csv.DictReader(fh_expected)  # Site	Class	Partition	Inserted?	Omega
            for expected_row in  expected_reader:
                site_1based = int(expected_row["Site"])
                expected_dnds = float(expected_row["Omega"])
                actual_dnds = actual_sites_dnds[site_1based-1]
                if round(expected_dnds, 2) != round(actual_dnds, 2):
                    print "Site: " + str(site_1based) + " ! expected=" + str(expected_dnds) + " actual=" + str(actual_dnds)
                else:
                    print "Site: " + str(site_1based) + " K  expected=" + str(expected_dnds) + " actual=" + str(actual_dnds)

        with open(expected_sites_dnds_csv, 'rU') as fh_expected:
            expected_reader = csv.DictReader(fh_expected)  # Site	Class	Partition	Inserted?	Omega
            for expected_row in  expected_reader:
                site_1based = int(expected_row["Site"])
                expected_dnds = float(expected_row["Omega"])
                actual_dnds = actual_sites_dnds[site_1based-1]


                self.assertAlmostEqual(expected_dnds, actual_dnds, 2,
                                       "Site: " + str(site_1based) + " expected=" + str(expected_dnds) + " actual=" + str(actual_dnds))

                print "Site: " + " OK " + str(site_1based) + " expected=" + str(expected_dnds) + " actual=" + str(actual_dnds)


if __name__ == "__main__":
    unittest.main()
