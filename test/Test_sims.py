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


    def test_prune(self):
        """
        Test edge cases for pruning
        :return:
        """
        treestr = "(((T1:0.1, T2:0.2)N3:0.35, (T3:0.3, T4:0.4)N4:0.45)N2:0.25, T5:0.5)N1;"
        tree = Phylo.read(StringIO(treestr), "newick")

        # Double check that the level order is correct
        actual_breadthfirst_nodes = [node.name for node in tree.find_clades(order="level")]
        expected_breadthfirst_nodes = ["N1", "N2", "T5", "N3", "N4", "T1", "T2", "T3", "T4"]
        self.assertListEqual(actual_breadthfirst_nodes, expected_breadthfirst_nodes)


        # TestCase:  parent is root, tip
        prune_idx =  expected_breadthfirst_nodes.index("T5")
        prune_clade = simulations.recombination.prune_by_idx(tree=tree, idx=prune_idx, order="level")

        self.assertEqual(prune_clade.name, "T5", msg="Expected T5 pruned but instead pruned " + str(prune_clade.name))
        expected_pruned_nodes = ["N2", "N3", "N4", "T1", "T2", "T3", "T4"]
        actual_pruned_nodes = [node.name for node in tree.find_clades(order="level")]
        self.assertListEqual(actual_pruned_nodes, expected_pruned_nodes,
                             msg="Actual=" + str(actual_pruned_nodes) +
                                 " Expected=" + str(expected_pruned_nodes))

        # TestCase:  parent is root, inner node
        tree = Phylo.read(StringIO(treestr), "newick")
        prune_idx =  expected_breadthfirst_nodes.index("N2")
        prune_clade = simulations.recombination.prune_by_idx(tree=tree, idx=prune_idx, order="level")

        self.assertEqual(prune_clade.name, "N2", msg="Expected N2 pruned but instead pruned " + str(prune_clade.name))
        expected_pruned_nodes = ["T5"]
        actual_pruned_nodes = [node.name for node in tree.find_clades(order="level")]
        self.assertListEqual(actual_pruned_nodes, expected_pruned_nodes,
                             msg="Actual=" + str(actual_pruned_nodes) +
                                 " Expected=" + str(expected_pruned_nodes))


        # TestCase:  parent is not root, inner node
        tree = Phylo.read(StringIO(treestr), "newick")
        prune_idx =  expected_breadthfirst_nodes.index("N3")
        prune_clade = simulations.recombination.prune_by_idx(tree=tree, idx=prune_idx, order="level")

        self.assertEqual(prune_clade.name, "N3", msg="Expected N2 pruned but instead pruned " + str(prune_clade.name))
        expected_pruned_nodes = ["N1", "N4", "T5", "T3", "T4"]
        actual_pruned_nodes = [node.name for node in tree.find_clades(order="level")]
        self.assertListEqual(actual_pruned_nodes, expected_pruned_nodes,
                             msg="Actual=" + str(actual_pruned_nodes) +
                                 " Expected=" + str(expected_pruned_nodes))

        # TestCase:  parent is not root, tip
        tree = Phylo.read(StringIO(treestr), "newick")
        prune_idx =  expected_breadthfirst_nodes.index("T4")
        prune_clade = simulations.recombination.prune_by_idx(tree=tree, idx=prune_idx, order="level")

        self.assertEqual(prune_clade.name, "T4", msg="Expected N2 pruned but instead pruned " + str(prune_clade.name))
        expected_pruned_nodes = ["N1", "N2", "T5", "N3", "T3", "T1", "T2"]
        actual_pruned_nodes = [node.name for node in tree.find_clades(order="level")]
        self.assertListEqual(actual_pruned_nodes, expected_pruned_nodes,
                             msg="Actual=" + str(actual_pruned_nodes) +
                                 " Expected=" + str(expected_pruned_nodes))




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


        # TestCase:  Check destination sister branch length where pruned clade collapsed destination sister parent
        treestr = "(((T1:0.1, T2:0.2)N3:0.35, (T3:0.3, T4:0.4)N4:0.45)N2:0.25, T5:0.5)N1;"
        tree = Phylo.read(StringIO(treestr), "newick")
        prune_idx = expected_breadthfirst_nodes.index("N3")
        dest_sister_clade = tree.find_any(name="N4")
        dest_new_par_br_len = 0.248459005204
        relocate_clade = simulations.recombination.prune_by_idx(tree=tree, idx=prune_idx, order="level")
        tree = simulations.recombination.graft(tree=tree, src_subtree=relocate_clade,
                                               dest_sister_clade=dest_sister_clade,
                                               dest_new_par_br_len=dest_new_par_br_len)

        exp_graft_str = "(((T3:0.3, T4:0.4)N4:0.451540995,(T1:0.1, T2:0.2)N3:0.35)N40:0.248459005204, T5:0.5)N1;"
        actual_strio = StringIO()
        Phylo.write(tree, actual_strio, "newick", format_branch_length='%1.9f')
        actual_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(exp_graft_str, actual_strio.getvalue())
        self.assertAlmostEqual(diff, 0, places=7,
                               msg="Grafted Tree not same as Expected Tree. Actual=" + actual_strio.getvalue() +
                                   " expected=" + exp_graft_str + " RF=" + str(diff))
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


    @staticmethod
    def get_parent(tree, target=None, **kwargs):
        """
        Helper method to get parent from tree and node
        :param Bio.Phylo.Tree tree:  tree
        :param Bio.Phylo.TreeMixin node:  node
        :return:
        """

        node_path = tree.get_path(target, **kwargs)

        if len(node_path) >= 2:
            parent = node_path[-2]
        else:
            parent = tree.root  #  parent is root

        return parent

    @staticmethod
    def get_sisters(tree, target=None, **kwargs):
        """
        Helper method to get sisters from tree and node
        :param Bio.Phylo.Tree tree:  tree
        :param Bio.Phylo.TreeMixin node:  node
        :return:
        """
        parent = TestSims.get_parent(tree=tree, target=target, **kwargs)
        curr_clade = tree.find_any(target, **kwargs)
        sisters = [clade for clade in parent.clades if clade != curr_clade]
        return sisters


    def test_recombo_tree(self):
        """
        Directly invoke recombination.recombo_tree() and check that the tips are the same but topology is different.
        Check that the grafted clade is in the right place.
        :return:
        """
        tmptree = tempfile.NamedTemporaryFile(mode="w+", delete=False)

        treestr = "(((T1:0.1, T2:0.2)N3:0.35, (T3:0.3, T4:0.4)N4:0.45)N2:0.25, T5:0.5)N1;"
        tmptree.write(treestr)
        tmptree.flush()  # Flush to buffer
        os.fsync(tmptree.file.fileno())  # flush to disk
        tmptree.close()


        # origtree = Phylo.read(StringIO(treestr), "newick")
        # expected_breadthfirst_nodes = ["N1", "N2", "T5", "N3", "N4", "T1", "T2", "T3", "T4"]
        # src_prune_id = expected_breadthfirst_nodes.index("T5")
        # dest_sister_clade = origtree.find_any(name="T3")
        # dest_new_par_br_len = 0.184629237455
        # # T5
        # # T3
        # # 0.184629237455
        #
        # graft_clade = simulations.recombination.prune_by_idx(tree=origtree, idx=src_prune_id, order="level")
        # recombo_tree = simulations.recombination.graft(tree=origtree,
        #                                                src_subtree=graft_clade,
        #                                                dest_sister_clade=dest_sister_clade,
        #                                                dest_new_par_br_len=dest_new_par_br_len)


        recombo_tree, graft_clade, dest_sister_clade, dest_new_par_br_len = simulations.recombination.make_recombo_tree(in_treefile=tmptree.name,
                                                            out_treefile=tmptree.name + ".recombo.nwk")


        print graft_clade.name
        print dest_sister_clade.name
        print dest_new_par_br_len

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
        Phylo.write(recombo_tree, recombo_strio, "newick", format_branch_length='%1.9f')
        recombo_strio.flush()
        diff = TestTopology.get_weighted_rf_dist_from_str(treestr, recombo_strio.getvalue())
        self.assertNotEqual(diff, 0, msg="Recombination Tree should be different from original tree. " +
                                      "Recombination Tree=" + recombo_strio.getvalue() +
                                      " Orig tree =" + treestr)


        # Check that grafted clade is grafted as a sister to dest_sister_clade
        recombo_dest_sis_parent = TestSims.get_parent(recombo_tree, target=graft_clade)
        sister_clade_names = [clade.name for clade in recombo_dest_sis_parent.clades]
        self.assertIn(dest_sister_clade.name, sister_clade_names,
                      msg="Expect clades " + str(dest_sister_clade.name) + " and " + str(graft_clade.name) + " to be sisters")


        # Check that destination sister clade's branch length has been updated
        orig_dest_sister = origtree.find_any(name=dest_sister_clade.name)
        # if the destination sister's original parent in the original tree was collapsed after pruning,
        # then the destination sister's branch is elongated by the original destination parent's branch length.
        orig_dest_sis_par = TestSims.get_parent(origtree, name=dest_sister_clade.name)
        if orig_dest_sis_par.branch_length is None:
            orig_dest_sis_par.branch_length = 0.0
        if recombo_tree.find_any(name=orig_dest_sis_par.name):
            expected_orig_sis_len = orig_dest_sister.branch_length - dest_new_par_br_len
        else:
            expected_orig_sis_len = orig_dest_sister.branch_length + orig_dest_sis_par.branch_length - dest_new_par_br_len

        self.assertEqual(expected_orig_sis_len, dest_sister_clade.branch_length,
                         msg="Expect destination sister clade to have branch length = " +
                             str(expected_orig_sis_len) +
                             " but got " + str(dest_sister_clade.branch_length) +
                             ".  Original destination sister branch len=" + str(orig_dest_sister.branch_length))


        # Check that robinson foulds distance is as expected
        expected_rf = 0.0

        orig_src_sister = TestSims.get_sisters(tree=origtree, name=graft_clade.name)[0]
        orig_src_par = TestSims.get_parent(tree=origtree, name=graft_clade.name)

        # Bipartition in original tree formed by source ancestors of pruned clade:
        # After pruning, the original source ancestors of the pruned clade will have fewer descendents
        # unless the pruned clade was regrafted as a descendent in the recombo tree
        orig_src_ancs = origtree.get_path(name=graft_clade.name)  # Traverse original ancestors of pruned clade, [just before root, curr clade]
        for orig_src_anc in orig_src_ancs[:-1]:
            # Check if the pruned clade was regrafted as another descendent in the recombo tree.
            # Add the ancestor clade's branch length if its descendents no longer include the pruned+regrafted clade
            # in the recombo tree.
            # if the pruned clade was regrafted as a sister in the recombo tree, handle it when we iterate through the destination sister ancestors.
            recombo_src_anc = recombo_tree.find_any(name=orig_src_anc.name)
            if recombo_src_anc and not recombo_src_anc.is_parent_of(name=graft_clade.name) and not orig_src_anc.name == dest_sister_clade.name:
                expected_rf += (orig_src_anc.branch_length or 0.0)
                expected_rf += (recombo_src_anc.branch_length or 0.0)


            # If original source parent only had 1 child after pruning, the original source parent collapses (ie no longer exists in recombo tree).
            # The original source sister will extend by the length of the parent in the recombo tree.
            # But the partition defined by the original source parent might still exist in the recombo tree if the pruned clade
            # was regrafted as a descendent of the original source sister.
            elif not recombo_src_anc and not orig_src_anc.is_parent_of(name=dest_sister_clade.name):
                expected_rf += 2 * (orig_src_anc.branch_length or 0.0)


        # Bipartition in recombination formed by ancestors of destination sister clade:
        # After regrafting, the original ancestors of the destination sister clade will have extra descendents
        # unless the regrafted clade was already a descendent.
        # Take on the root of the recombo tree since it's possible that the pruned clade is a child of the root, which causes
        # the original root to collapse and the original source sister to become the new root in the recombo tree
        recombo_dest_ancs = [recombo_tree.root] + recombo_tree.get_path(name=dest_sister_clade.name)   # Traverse original ancestors of destination sister clade
        for recombo_dest_anc in recombo_dest_ancs[:-1]:
            # We've already accounted bipartitions formed by the ancestors of the grafted clade
            orig_dest_anc = origtree.find_any(name=recombo_dest_anc.name)
            if orig_dest_anc and not orig_dest_anc.is_parent_of(name=graft_clade.name):
                # If the recombo destination ancestor is the original source sister, then we need to take into account
                # that the the original source parent was collapsed and its branch length added onto the original source sister.
                # Also the bipartition formed by original source parent is now formed by by the original sister.
                if orig_dest_anc.name == orig_src_sister.name:
                    expected_rf += (orig_dest_anc.branch_length or 0.0)
                    expected_rf += abs((orig_src_par.branch_length or 0.0) - (recombo_dest_anc.branch_length or 0.0))
                else:
                    # Add the original destination ancestor clade's branch length if its descendents
                    # didn't include the pruned+regrafted clade in the original tree.
                    expected_rf += (orig_dest_anc.branch_length or 0.0)
                    expected_rf += (recombo_dest_anc.branch_length or 0.0)
            elif not orig_dest_anc:
                # New intermediate node as new parent of recombo tree destination sister clade
                # If the destination sister is an ancestor of the source,
                # then the bipartition formed by the destination sister is now formed by the new intermediate node.
                if orig_dest_sister.is_parent_of(name=graft_clade.name):
                    expected_rf += 2*abs((orig_dest_sister.branch_length or 0.0) - (recombo_dest_anc.branch_length or 0.0))
                else:
                    expected_rf += 2*(recombo_dest_anc.branch_length or 0.0)


        self.assertAlmostEqual(expected_rf, diff, places=7, msg="Expect rf=" + str(expected_rf) + " but got " + str(diff) +
                                                                " RecomboTree=" + recombo_strio.getvalue())


        # cleanup
        recombo_strio.close()
        if os.path.exists(tmptree.name):
            os.remove(tmptree.name)
        if os.path.exists(tmptree.name + ".recombo.nwk"):
            os.remove(tmptree.name + ".recombo.nwk")


    @staticmethod
    def hamming_distance(s1, s2):
        """
        Gets number of bases different between sequences
        :return int:
        """
        diffs = 0
        for i, letter1 in enumerate(s1):
            if i >= len(s2):
                diffs += 1
            elif letter1 != s2[i]:
                diffs += 1

        return diffs



    def test_reconstruct_full_popn(self):
        """
        How faithfully does HyPhy reconstruct the INDELible ancestral sequences given
        the topology and FastTree reconstructed branch lengths, and INDELible tip sequences?
        :return:
        """

        indelible_anc_fasta = SIM_DATA_DIR + os.sep + "fullpopn" + os.sep + SIM_DATA_FILENAME_PREFIX + "_ANCESTRAL.fasta"
        indelible_aln =SeqIO.to_dict(SeqIO.parse(indelible_anc_fasta, "fasta"))


        # If there are recombination breaks, just take the first section
        hyphy_anc_fasta = glob.glob(SIM_DATA_DIR + os.sep + "subs" + os.sep + SIM_DATA_FILENAME_PREFIX + ".break.1_*.anc.fasta")[0]
        hyphy_aln = SeqIO.to_dict(SeqIO.parse(hyphy_anc_fasta, "fasta"))

        # compare hyphy sequence against indelible sequence
        total_dist = 0.0
        total_seq = 0.0

        for anc_id in hyphy_aln.keys():
            if anc_id.find("otu") >= 0:
                continue  # Don't include tips in calculations
            hyphy_seq = hyphy_aln[anc_id].seq
            indelible_seq = indelible_aln[anc_id].seq
            dist = TestSims.hamming_distance(hyphy_seq, indelible_seq)
            total_dist += dist
            total_seq += 1

            print "Ancestor ID = " + anc_id + " dist=" + str(dist)

        print "average hamming distance = " + str(total_dist) + "/" + str(total_seq) + " = " + str(total_dist/total_seq)




if __name__ == "__main__":
    unittest.main()
