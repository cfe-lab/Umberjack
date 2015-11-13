# Tests that the simulated data are actually simulating what they're supposed to

# calculates the subsitution rate at the nucleotide level and amino acid level
# checks that the expected substitution rate is the actual substitution rate
import unittest
from Bio import Phylo
from Bio import Seq as Seq
import simulations.indelible.indelible_handler as indelibler
from cStringIO import StringIO
import csv
import os
import subprocess
import config.settings
import shutil
import ConfigParser
import Bio.SeqIO as SeqIO
from Bio.Align import MultipleSeqAlignment
import simulations.recombination
import tempfile
import glob
from cStringIO import StringIO
from test_topology import TestTopology

# Simulation Configs
SIM_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations"
SIM_BIN_DIR = SIM_DIR + os.sep + "bin"
SIM_DATA_FILENAME_PREFIX = "umberjack_unittest"
SIM_DATA_DIR = SIM_DIR + os.sep + "data" + os.sep + SIM_DATA_FILENAME_PREFIX
SIM_DATA_CONFIG_FILE = SIM_DATA_DIR + os.sep + "umberjack_unittest_sim.conf"

# Executables
SIM_PIPELINE_PY = os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "sim_pipeline.py"

# Config File section
SECTION = "sim"

class TestSims(unittest.TestCase):

    def setUp(self):
        """
        Generate simulated data for unit tests.
        The sim_pipeline.py should be smart enough not to regenerate files if they already exist to save time.
        """
        config.settings.setup_logging()
        self.configparser = ConfigParser.RawConfigParser()
        self.configparser.read(SIM_DATA_CONFIG_FILE)
        subprocess.check_call(["python", SIM_PIPELINE_PY, SIM_DATA_CONFIG_FILE])


    @staticmethod
    def get_tree_len(treefilename):
        """
        Gets the tree length, i.e.  the sum of every branch in the tree.  Excludes root branch length.
        """
        tree = Phylo.read(treefilename, "newick")
        root_branch_length = tree.clade.branch_length if tree.clade.branch_length is not None else 0
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
        """
        Check that INDELible trees lengths are scaled to the rate given in config file.
        :return:
        """

        # for each scaling rate, check that indelible scaled the tree correctly
        for expected_tree_len in [int(x) for x in self.configparser.get(SECTION, "INDELIBLE_SCALING_RATES").split(",")]:

            indelible_tree_txt = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.txt"
            indelible_tree_str = indelibler.get_tree_string(indelible_tree_txt)

            indelible_tree_handle = StringIO(indelible_tree_str)
            actual_tree_len = TestSims.get_tree_len(indelible_tree_handle)
            indelible_tree_handle.close()

            self.assertAlmostEqual(expected_tree_len, actual_tree_len, 0,
                               "Expected tree len=" + str(expected_tree_len) + " but got " + str(actual_tree_len))


    def test_sample_genomes(self):
        """
        Tests that the final population sequences are comprised of the correct blocks of
        INDELible population sequences at different scaling rates.
        :return:
        """
        # data/umberjack_unittest/mixed/umberjack_unittest.mixed.fasta
        filename_prefix = self.configparser.get(SECTION, "FILENAME_PREFIX")
        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + filename_prefix + ".mixed.fasta"
        final_popn_rates_csv = SIM_DATA_DIR + os.sep + "mixed" + os.sep + filename_prefix + ".mixed.rates.csv"



        # Each scaling factor should be represented equally
        site_to_scale = []
        scale_to_codons = dict()
        with open(final_popn_rates_csv, "rU") as fh_in:
            # Site	Scaling_factor	Rate_class	Omega
            reader = csv.DictReader(fh_in)
            for row in reader:
                codonsite_0based = int(row["Site"]) - 1
                scaling_factor = row["Scaling_factor"]
                site_to_scale.append((codonsite_0based, scaling_factor))
                if scaling_factor not in  scale_to_codons:
                    scale_to_codons[scaling_factor] = [codonsite_0based]
                else:
                    scale_to_codons[scaling_factor].extend([codonsite_0based])

        codons_perblock = self.configparser.getint(SECTION, "CODONS_PER_BLOCK")
        total_codons = self.configparser.getint(SECTION, "NUM_CODON_SITES")
        scaling_factors = [int(x) for x in self.configparser.get(SECTION, "INDELIBLE_SCALING_RATES").split(",")]
        expected_codons_per_scale = total_codons/len(scaling_factors)

        actual_codon_counts_per_scale = [len(x) for x in scale_to_codons.values()]
        self.assertItemsEqual([expected_codons_per_scale] * len(scaling_factors),
                              actual_codon_counts_per_scale,
                              "Expect same number of codons per scaling factor")


        # Each scaling factor should be in blocks of consecutive codons
        site_to_scale = sorted(site_to_scale, key=lambda (site, scale): site)
        last_scale = None
        for codonsite_0based, scale in site_to_scale:
            # We are not at the beginning of a new block
            if codonsite_0based % codons_perblock != 0 and last_scale is not None:
                # if this site's scaling factor isn't the same as the last scaling factor, crap out
                self.assertEqual(last_scale, scale, "Expect that consecutive codons in block should have same scaling factor")
            last_scale = scale


        # Each block at a certain scaling rate in the final population should be equal
        # to the same block in the original INDELibel population at the same scaling rate.
        # Check for both tip fastas and ancestral fastas.
        actual_aln = MultipleSeqAlignment(SeqIO.parse(final_popn_fasta, "fasta"))
        for scale, codons in scale_to_codons.iteritems():
            # data/umberjack_unittest/50/umberjack_unittest.50.fasta
            orig_popn_fasta = SIM_DATA_DIR + os.sep + str(scale) + os.sep + filename_prefix + ".{}.fasta".format(scale)
            expected_aln = MultipleSeqAlignment(SeqIO.parse(orig_popn_fasta, "fasta"))
            for codonsite_0based in codons:
                nuc_start_0based = codonsite_0based * 3
                for nuc_0based in range(nuc_start_0based, nuc_start_0based + 3):
                    self.assertEqual(expected_aln[:, nuc_0based],
                                     actual_aln[:, nuc_0based],
                                     "Expect nucleotide site (0based) " + str(nuc_0based) + " in " + final_popn_fasta +
                                     " to originate from " + orig_popn_fasta)

        # inner node (ancestral sequences)
        # Be careful!  HyPhy batch scripts will auto write its reconstructed ancestors and tips to .mixed.anc.fasta,
        # which is not the same as the concatenation of original ancestors, scaled at different rates in .mixed.ancestral.fasta
        final_anc_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + filename_prefix + ".mixed.ancestral.fasta"
        actual_anc_aln = MultipleSeqAlignment(SeqIO.parse(final_anc_popn_fasta, "fasta"))

        for scale, codons in scale_to_codons.iteritems():

            orig_anc_popn_fasta = SIM_DATA_DIR + os.sep + str(scale) + os.sep + filename_prefix + ".{}_ANCESTRAL.fasta".format(scale)
            expected_anc_aln = MultipleSeqAlignment(SeqIO.parse(orig_anc_popn_fasta, "fasta"))

            for codonsite_0based in codons:
                nuc_start_0based = codonsite_0based * 3
                for nuc_0based in range(nuc_start_0based, nuc_start_0based + 3):
                    self.assertEqual(expected_anc_aln[:, nuc_0based],
                                     actual_anc_aln[:, nuc_0based],
                                     "Expect nucleotide site (0based) " + str(nuc_0based) + " in " + final_anc_popn_fasta +
                                     " to originate from " + orig_anc_popn_fasta)



    def test_prune_regraft(self):
        """
        Tests that the graft inserts grafted clade in correct position
        :return:
        """
        treestr = "(((T1:0.1, T2:0.2)N3:0.35, (T3:0.3, T4:0.4)N4:0.45)N2:0.25, T5:0.5)N1;"
        tree = Phylo.read(StringIO(treestr), "newick")

        # Double check that the level order is correct
        actual_breadthfirst_nodes = [node.name for node in tree.find_clades(order="level")]
        expected_breadthfirst_nodes = ["N1", "N2", "T5", "N3", "N4", "T1", "T2", "T3", "T4"]
        self.assertListEqual(actual_breadthfirst_nodes, expected_breadthfirst_nodes)

        # Check that the prune_by_idx() actually pruned the right node
        relocate_idx =  expected_breadthfirst_nodes.index("N3")
        graft_clade = simulations.recombination.prune_by_idx(tree=tree, idx=relocate_idx, order="level")

        # Expect that pruning clade N3 causes N2 to be a singleton node.
        # N2 is then collapsed so that N1 children are N4 and T5
        expected_prune_str = "((T3:0.3, T4:0.4)N4:0.70, T5:0.5)N1;"
        actual_strio = StringIO()
        Phylo.write(tree, actual_strio, "newick")
        actual_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(expected_prune_str, actual_strio.getvalue())
        self.assertEqual(diff, 0, msg="Pruned Tree not same as Expected Tree. Actual=" + actual_strio.getvalue() +
                                             " expected=" + expected_prune_str)
        actual_strio.close()

        # Check that regrafting the pruned clade inserts it in the correct location
        dest_sister_clade = tree.find_clades(name="T5").next()
        dest_new_par_br_len = 0.2
        tree = simulations.recombination.graft(tree=tree, src_subtree=graft_clade,
                                               dest_sister_clade=dest_sister_clade,
                                               dest_new_par_br_len=dest_new_par_br_len)

        exp_graft_str = "((T3:0.3, T4:0.4)N4:0.70, ((T1:0.1, T2:0.2)N3:0.35, T5:0.3):0.2)N1;"
        actual_strio = StringIO()
        Phylo.write(tree, actual_strio, "newick")
        actual_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(exp_graft_str, actual_strio.getvalue())

        self.assertEqual(diff, 0, msg="Grafted Tree not same as Expected Tree. Actual=" + actual_strio.getvalue() +
                                             " expected=" + exp_graft_str)

        actual_strio.close()



    def test_prune_regraft_polytomy(self):
        """
        Check that regrafting the pruned clade on top of existing node
        inserts the pruned clade as a polytomy sister of the existing node
        :return:
        """
        treestr = "(((T1:0.1, T2:0.2)N3:0.35, (T3:0.3, T4:0.4)N4:0.45)N2:0.25, T5:0.5)N1;"
        tree = Phylo.read(StringIO(treestr), "newick")

        # Double check that the level order is correct
        expected_breadthfirst_nodes = ["N1", "N2", "T5", "N3", "N4", "T1", "T2", "T3", "T4"]

        # Check that the prune_by_idx() actually pruned the right node
        relocate_idx =  expected_breadthfirst_nodes.index("N3")
        graft_clade = simulations.recombination.prune_by_idx(tree=tree, idx=relocate_idx, order="level")

        # Expect that pruning clade N3 causes N2 to be a singleton node.
        # N2 is then collapsed so that N1 children are N4 and T5
        expected_prune_str = "((T3:0.3, T4:0.4)N4:0.70, T5:0.5)N1;"
        actual_strio = StringIO()
        Phylo.write(tree, actual_strio, "newick")
        actual_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(expected_prune_str, actual_strio.getvalue())
        self.assertEqual(diff, 0, msg="Pruned Tree not same as Expected Tree. Actual=" + actual_strio.getvalue() +
                                             " expected=" + expected_prune_str)
        actual_strio.close()

        # Check that regrafting the pruned clade inserts it in the correct location
        dest_sister_clade = tree.find_clades(name="T5").next()
        dest_new_par_br_len = 0.0
        tree = simulations.recombination.graft(tree=tree, src_subtree=graft_clade,
                                               dest_sister_clade=dest_sister_clade,
                                               dest_new_par_br_len=dest_new_par_br_len)

        exp_graft_str = "((T3:0.3, T4:0.4)N4:0.70, (T1:0.1, T2:0.2)N3:0.35, T5:0.5)N1;"
        actual_strio = StringIO()
        Phylo.write(tree, actual_strio, "newick")
        actual_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(exp_graft_str, actual_strio.getvalue())

        self.assertEqual(diff, 0, msg="Grafted Tree not same as Expected Tree. Actual=" + actual_strio.getvalue() +
                                             " expected=" + exp_graft_str)

        actual_strio.close()


    def test_regraft_root(self):
        """
        Check that you can't graft a clade as the new root
        :return:
        """

        treestr = "((T3:0.3, T4:0.4)N4:0.70, T5:0.5)N1;"
        tree = Phylo.read(StringIO(treestr), "newick")
        dest_sister_clade = tree.find_clades(name="T5").next()
        dest_new_par_br_len = 0.2

        subtree_str = "(T3:0.3, T4:0.4)N4:0.45;"
        subtree = Phylo.read(StringIO(subtree_str), "newick")

        self.assertRaises(ValueError, simulations.recombination.graft,
                          **dict(tree=tree, src_subtree=subtree,
                                 dest_sister_clade=dest_sister_clade, dest_new_par_br_len=dest_new_par_br_len))



    def test_recombination(self):
        """
        Test that simulating a recombination tree by randomly pruning and regraphing a clade works.
        :return:
        """

        prev_brlen = None
        prev_tips = None
        prev_treefiles = []
        for recombo_treefile in glob.glob(SIM_DATA_DIR + os.sep + "topology" + os.sep + "*.break.*.rename.nwk"):
            tree = Phylo.read(recombo_treefile, "newick")

            # Check that total branch length in each tree is the same in the umberjack unittest recombination trees
            curr_brlen = tree.total_branch_length()
            if prev_brlen is not None:
                self.assertAlmostEqual(prev_brlen, curr_brlen,
                                       msg="Expect all recombination trees to have same total branch length.  " +
                                           "Difference in {}={} and {}={}".format(recombo_treefile, curr_brlen,
                                                                                  prev_treefiles[-1], prev_brlen))

            # Check that all the same tips exist in each tree
            curr_tips = [tip.name for tip in tree.get_terminals()]
            if prev_tips is not None:
                self.assertItemsEqual(prev_tips, curr_tips,
                                      msg="Expect all recombination trees to have same tips.  " +
                                           "Difference in {}={} and {}={}".format(recombo_treefile, str(prev_tips),
                                                                                  prev_treefiles[-1], str(curr_tips)))


            # Check that each recombination tree is different from all other treefiles
            for prev_treefile in prev_treefiles:
                diff = TestTopology.get_weighted_rf_dist(recombo_treefile, prev_treefile)
                self.assertNotEqual(diff, 0, msg="Recombination Tree is same as another tree. Current Tree=" + recombo_treefile +
                                                 " Previous Tree =" + prev_treefile)

            prev_brlen  = curr_brlen
            prev_tips = curr_tips
            prev_treefiles.append(recombo_treefile)


    def test_recombo_tree(self):
        """
        Directly invoke recombination.recombo_tree() and check that the tips are the same but topology is different.
        :return:
        """
        tmptree = tempfile.NamedTemporaryFile(mode="w+", delete=False)

        treestr = "(((T1:0.1, T2:0.2)N3:0.35, (T3:0.3, T4:0.4)N4:0.45)N2:0.25, T5:0.5)N1;"
        tmptree.write(treestr)
        tmptree.flush()  # Flush to buffer
        os.fsync(tmptree.file.fileno())  # flush to disk
        tmptree.close()

        graft_clade, recombo_tree = simulations.recombination.make_recombo_tree(in_treefile=tmptree.name,
                                                             out_treefile=tmptree.name + ".recombo.nwk",
                                                             seed=20)

        # Check that the tree tips are the same
        origtree = Phylo.read(StringIO(treestr), "newick")
        orig_tips = [tip.name for tip in origtree.get_terminals()]
        recombo_tips = [tip.name for tip in recombo_tree.get_terminals()]
        self.assertItemsEqual(orig_tips, recombo_tips,
                              msg="Recombination Tree should have same tips as the original tree.  " +
                                  "Recombo tips=" + str(recombo_tips) +
                                  "Orig tree tips=" + str(orig_tips))

        # Check that the topology is different
        recombo_strio = StringIO()
        Phylo.write(recombo_tree, recombo_strio, "newick")
        recombo_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(treestr, recombo_strio.getvalue())
        self.assertNotEqual(diff, 0, msg="Recombination Tree should be different from original tree. " +
                                      "Recombination Tree=" + recombo_strio.getvalue() +
                                      " Orig tree =" + treestr)
        recombo_strio.close()

        # cleanup
        if os.path.exists(tmptree.name):
            os.remove(tmptree.name)
        if os.path.exists(tmptree.name + ".recombo.nwk"):
            os.remove(tmptree.name + ".recombo.nwk")





if __name__ == "__main__":
    unittest.main()
