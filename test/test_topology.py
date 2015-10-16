"""
Checks that topology is stable from FastTree across genome.
- Test that FastTree reconstrucsts tree with same topology as coalescent simulated topology using INDELible fasta sequences.
- Test that INDELible tree topology matches the topology of simulated coalescent
- Test that FastTree reconstruction of tree from subsampling without replacement from the final full population
    matches coalescent topology with unsampled tips pruned
"""
import unittest
from Bio import Phylo
from Bio import Seq as Seq
import simulations.indelible.indelible_handler as indelibler
import simulations.phylo_dnds as phylo_dnds
from cStringIO import StringIO
import csv
import os
import subprocess
import config.settings
import shutil
import ConfigParser
import Bio.SeqIO as SeqIO
from Bio.Align import MultipleSeqAlignment
import fasttree.fasttree_handler
import random
from collections import defaultdict

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

class TestTopology(unittest.TestCase):

    def setUp(self):
        """
        Generate simulated data for unit tests unless it already exists (takes a long time to run)
        """
        config.settings.setup_logging()

        self.configparser = ConfigParser.RawConfigParser()
        self.configparser.read(SIM_DATA_CONFIG_FILE)
        filename_prefix = self.configparser.get(SECTION, "FILENAME_PREFIX")


        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + filename_prefix + ".mixed.fasta"
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#
        if not os.path.exists(final_popn_fasta) or not os.path.getsize(final_popn_fasta) or not os.path.exists(final_popn_tree_file):
            subprocess.check_call(["python", SIM_PIPELINE_PY, SIM_DATA_CONFIG_FILE])


        # Create newick trees of the final full population and the original coalescent without branch lengths to compare in splitstree
        # /data/umberjack_unittest/mixed/umberjack_unittest.mixed.nwk
        final_popn_tree_nolen_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nolen.nwk"#


        final_popn_tree = Phylo.read(final_popn_tree_file, "newick")
        for clade in final_popn_tree.find_clades():
            clade.branch_length = 1
        Phylo.write(final_popn_tree, final_popn_tree_nolen_file, "newick")

        # /data/umberjack_unittest/umberjack_unittest.rename.nwk
        orig_tree_file = SIM_DATA_DIR + os.sep + SIM_DATA_FILENAME_PREFIX + ".rename.nwk"
        orig_tree_nolen_file = SIM_DATA_DIR + os.sep + SIM_DATA_FILENAME_PREFIX + ".rename.nolen.nwk"
        orig_tree = Phylo.read(orig_tree_file, "newick")
        for clade in orig_tree.find_clades():
            clade.branch_length = 1
        Phylo.write(orig_tree, orig_tree_nolen_file, "newick")

        try:
            import rpy2.robjects
        except ImportError:
            print "Unable to import rpy2 and thus unable to check tree toplogies"
            raise



    def test_full_popn_tree(self):
        """
        Tests that the FastTree tree of the full population uses the same topology as the input tree generated from the coalescent simulator.
        :return:
        """
        # Use R phangorn package to calculate
        import rpy2.robjects as ro

        # /data/umberjack_unittest/mixed/umberjack_unittest.mixed.nwk
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#

        # /data/umberjack_unittest/umberjack_unittest.rename.nwk
        orig_tree_file = SIM_DATA_DIR + os.sep + SIM_DATA_FILENAME_PREFIX + ".rename.nwk"

        ro.r("library(phangorn)")
        ro.r("library(ape)")
        ro.r("orig_tree <- read.tree('{}')".format(orig_tree_file))
        ro.r("reconstruct_tree <- read.tree('{}')".format(final_popn_tree_file))
        # robinson foulds distance
        ro.r("rf_dist <- RF.dist(orig_tree, reconstruct_tree)")
        rf_dist = ro.r("rf_dist")[0]


        self.assertEqual(rf_dist, 0,
                         "Final Full Population Tree should have same topology as original coalescent tree, " +
                         " Instead robinson foulds distance =" + str(rf_dist))



    def test_indelible_topology(self):
        """
        Test that the INDELible scaled tree has the same topology as the original coalescent tree
        :return:
        """
        import rpy2.robjects as ro

        for expected_tree_len in [int(x) for x in self.configparser.get(SECTION, "INDELIBLE_SCALING_RATES").split(",")]:
            # INDELible creates a TSV with the tree.  Convert to separate newick file.
            indelible_tree_txt = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.txt"
            indelible_tree_io = indelibler.get_tree_stringio(indelible_tree_txt)
            indelible_tree = Phylo.read(indelible_tree_io, "newick")
            indelible_nwk = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.nwk"
            Phylo.write(indelible_tree, indelible_nwk, "newick")
            indelible_tree_io.close()


            # /data/umberjack_unittest/umberjack_unittest.rename.nwk
            orig_tree_file = SIM_DATA_DIR + os.sep + SIM_DATA_FILENAME_PREFIX + ".rename.nwk"


            ro.r("library(phangorn)")
            ro.r("library(ape)")
            ro.r("orig_tree <- read.tree('{}')".format(orig_tree_file))
            ro.r("reconstruct_tree <- read.tree('{}')".format(indelible_nwk))
            # robinson foulds distance
            ro.r("rf_dist <- RF.dist(orig_tree, reconstruct_tree)")
            rf_dist = ro.r("rf_dist")[0]


            self.assertEqual(rf_dist, 0,
                         "INDELible  tree " + indelible_nwk + " should have same topology as coalescent tree " +
                         orig_tree_file + ". Instead robinson foulds distance =" + str(rf_dist))



    def test_fasttree_topology(self):
        """
        Test that FastTree can reproduce the same topology of an INDELible tree given the INDELibel fasta sequences as input
        and the INDELIble tree as a constraint.
        We do this on the full genome population scaled at only 1 mutation scaling rate
        so that we don't have multiple mutation rates per sequence confounding issues with tree reconstruction.
        :return:
        """

        import rpy2.robjects as ro

        for expected_tree_len in [int(x) for x in self.configparser.get(SECTION, "INDELIBLE_SCALING_RATES").split(",")]:
            # INDELible creates a TSV with the tree.  Convert to separate newick file.
            indelible_tree_txt = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.txt"
            indelible_tree_io = indelibler.get_tree_stringio(indelible_tree_txt)
            indelible_tree = Phylo.read(indelible_tree_io, "newick")
            indelible_nwk = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.nwk"
            Phylo.write(indelible_tree, indelible_nwk, "newick")
            indelible_tree_io.close()


            # Try to reproduce the INDELible tree with the INDELible fasta sequences as input and the INDELible tree as topology constraint
            indelible_tip_fasta = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "{}.{}_TRUE.fasta".format(SIM_DATA_FILENAME_PREFIX, expected_tree_len)
            repro_treefile = fasttree.fasttree_handler.make_tree_repro(fasta_fname=indelible_tip_fasta,
                                                                       intree_fname=indelible_nwk)


            ro.r("library(phangorn)")
            ro.r("library(ape)")
            ro.r("orig_tree <- read.tree('{}')".format(indelible_nwk))
            ro.r("reconstruct_tree <- read.tree('{}')".format(repro_treefile))
            # robinson foulds distance
            ro.r("rf_dist <- RF.dist(orig_tree, reconstruct_tree)")
            rf_dist = ro.r("rf_dist")[0]


            self.assertEqual(rf_dist, 0,
                         "Reconstructed FastTree tree " + repro_treefile + " should have same topology as input tree " +
                         indelible_nwk + ". Instead robinson foulds distance =" + str(rf_dist))


    @staticmethod
    def subsample(treefile, fastafile, out_treefile, out_fastafile, fraction, replace=False, seed=None):
        """
        Subsamples population with or without replacement up to the desired fraction.  Returns subsampled tree and fasta.
        :param str treefile:  filepath to newick population tree
        :param str fastafile:  filepath to population sequences in fasta format
        :param float fraction:  fraction of population to sample
        :param bool replace:  whether to sample with replacement
        :param int seed:  the random seed to randomly select individuals.
        :return [str]:  list of sequence names selected.  If multiple copies of same sequence, then rename "<name>_<copy">
        """
        tree = Phylo.read(treefile, "newick")
        tips = tree.get_terminals()

        sample_size = int(round(fraction *  len(tips)))

        randomizer = random.Random(seed)
        if not replace:
            selecttip_to_count = dict((tip.name, 1) for tip in randomizer.sample(tips, sample_size))
        else:
            selecttip_to_count = defaultdict(int)
            for i in range(0, sample_size):
                tip = randomizer.choice(tips)
                selecttip_to_count[tip.name] += 1

        tipname_sample = []
        seqdict = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
        with open(out_fastafile, 'w') as fh_out_fasta:
            for tipname, count in selecttip_to_count.iteritems():
                if count >= 2:
                    for i in range(0, count):
                        new_tipname = "{}_{}".format(tipname, i)
                        fh_out_fasta.write(">" + new_tipname + "\n")
                        fh_out_fasta.write(str(seqdict[tipname].seq) + "\n")
                        tipname_sample.extend([new_tipname])
                elif count == 1:
                    fh_out_fasta.write(">{}\n".format(tipname))
                    fh_out_fasta.write(str(seqdict[tipname].seq) + "\n")
                    tipname_sample.extend([tipname])

        # Prune the tree
        for tip in tips:
            if tip.name not in selecttip_to_count.keys() or selecttip_to_count[tip.name] == 0:
                tree.prune(tip)

        # Rename tips that have been sampled multiple times to the name of the first copy.
        for tipname, count in selecttip_to_count.iteritems():
            if count >= 2:
                tip = tree.find_clades(target=tipname, terminal=True).next()
                tip.name = "{}_0".format(tipname)

        Phylo.write(tree, out_treefile, "newick")

        return tipname_sample


    def test_subsample(self):
        """
        Test that FastTree reconstruction of subsampled INDELible population has same topology as INDELible tree.
        Don't use final population whose sequences are concatenated from multiple muitation rates.
        Instead, just use a population whose entire genome has been scaled by a single mutation rate.

        """
        import rpy2.robjects as ro

        seed = self.configparser.getint(SECTION, "SEED")
        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        for expected_tree_len in [int(x) for x in self.configparser.get(SECTION, "INDELIBLE_SCALING_RATES").split(",")]:
            for sample_fraction in [0.5]:
                # INDELible creates a TSV with the tree.  Convert to separate newick file.
                indelible_tree_txt = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.txt"
                indelible_tree_io = indelibler.get_tree_stringio(indelible_tree_txt)
                indelible_tree = Phylo.read(indelible_tree_io, "newick")
                indelible_nwk = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.nwk"
                Phylo.write(indelible_tree, indelible_nwk, "newick")
                indelible_tree_io.close()


                # Make an expected tree by subsampling the indelibel tree
                # Subsample the indelible population sequences and feed that as input to FastTree for reconstruction
                indelible_tip_fasta = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "{}.{}_TRUE.fasta".format(SIM_DATA_FILENAME_PREFIX, expected_tree_len)
                subsample_fasta = indelible_tip_fasta.replace(".fasta", ".prune.{:.1}.fasta".format(sample_fraction))
                expected_subsample_treefile = indelible_nwk.replace(".nwk", ".prune.{:.1}.nwk".format(sample_fraction))

                TestTopology.subsample(treefile=indelible_nwk, fastafile=indelible_tip_fasta,
                                       out_treefile=expected_subsample_treefile, out_fastafile=subsample_fasta,
                                       seed=seed, fraction=sample_fraction)



                # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
                repro_treefile = fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                                     threads=fasttree_threads, debug=True)

                ro.r("library(phangorn)")
                ro.r("library(ape)")
                ro.r("expected_tree <- read.tree('{}')".format(expected_subsample_treefile))
                ro.r("reconstruct_tree <- read.tree('{}')".format(repro_treefile))

                # FastTree makes unrooted tree; its root is arbitrarily set.  Reroot trees with same root before compare.
                # Use first sequence in subsampled fasta as the root for both trees.
                first = SeqIO.parse(subsample_fasta, "fasta").next()

                ro.r("reconstruct_tree_rooted <- root(reconstruct_tree, outgroup='{}', resolve.root=TRUE)".format(first.id))
                ro.r("expected_tree_rooted <- root(expected_tree, outgroup='{}', resolve.root=TRUE)".format(first.id))

                # robinson foulds distance
                ro.r("rf_dist <- RF.dist(expected_tree_rooted, reconstruct_tree_rooted)")
                rf_dist = ro.r("rf_dist")[0]


                self.assertEqual(rf_dist, 0,
                             "Reconstructed FastTree tree " + repro_treefile + " should have same topology as input tree " +
                             expected_subsample_treefile + ". Instead robinson foulds distance =" + str(rf_dist))


    @staticmethod
    def prune_copies(treefile, out_treefile):
        """
        Finds polytomies made from copies of the same sequence indicated by tip names with format "[name]_[copy number]".
        Prunes all but the first copy "[name]_0" when it encounters these copies.
        If a tip is in a separate clade as the rest of the copies, it doesn't not prune it.
        :return:
        """
        repro_tree = Phylo.read(treefile, 'newick')
        for tip in repro_tree.get_terminals():
            copy_index = tip.name.find("_0")
            tipname_base = tip.name[:copy_index+1]
            if copy_index >= 0:  # multiple copies if suffixed with "_<copy number>"
                parent = repro_tree.get_path(tip)[-2]
                sisters = list(parent.clades)  # make a shallow copy of the sisters so when we prune them, we still iterate through all remaining sisters
                for sister in sisters:
                    if sister.name == tip.name:
                        continue
                    elif sister.name.find(tipname_base) >= 0:  # another copy of the same tip
                        repro_tree.prune(sister)

        Phylo.write(repro_tree, out_treefile, "newick")


    def test_subsample_replace(self):
        """
        Test that FastTree reconstruction of subsampled with replacement INDELible population has same topology as INDELible tree
        that has been pruned for unsampled tips.
        Don't use final population whose sequences are concatenated from multiple muitation rates.
        Instead, just use a population whose entire genome has been scaled by a single mutation rate.
        :return:
        """
        import rpy2.robjects as ro

        seed = self.configparser.getint(SECTION, "SEED")
        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        for expected_tree_len in [int(x) for x in self.configparser.get(SECTION, "INDELIBLE_SCALING_RATES").split(",")]:
            for sample_fraction in [0.5]:
                # INDELible creates a TSV with the tree.  Convert to separate newick file.
                indelible_tree_txt = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.txt"
                indelible_tree_io = indelibler.get_tree_stringio(indelible_tree_txt)
                indelible_tree = Phylo.read(indelible_tree_io, "newick")
                indelible_nwk = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "trees.nwk"
                Phylo.write(indelible_tree, indelible_nwk, "newick")
                indelible_tree_io.close()


                # Make an expected tree by subsampling the indelible tree.
                # Subsample the indelible population sequences and feed that as input to FastTree for reconstruction.
                # Since we allow replacement, add copies of the sequence with suffix "_<copy number>" to the fasta when necessary.
                indelible_tip_fasta = SIM_DATA_DIR + os.sep + str(expected_tree_len) + os.sep + "{}.{}_TRUE.fasta".format(SIM_DATA_FILENAME_PREFIX, expected_tree_len)
                subsample_fasta = indelible_tip_fasta.replace(".fasta", ".resample.{:.1}.fasta".format(sample_fraction))
                expected_subsample_treefile = indelible_nwk.replace(".nwk", ".resample.{:.1}.nwk".format(sample_fraction))

                TestTopology.subsample(treefile=indelible_nwk, fastafile=indelible_tip_fasta,
                                       out_treefile=expected_subsample_treefile, out_fastafile=subsample_fasta,
                                       seed=seed, fraction=sample_fraction, replace=True)



                # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
                repro_treefile = fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                                     threads=fasttree_threads, debug=True)

                # robinson foulds distance only well defined for binary trees.  Prune sister copies of the same sequence from tree.
                # If the copies of the same sequence aren't polytomies in the tree, then FastTree screwed up.

                repro_binary_treefile =  repro_treefile.replace(".nwk", ".binary.nwk")
                TestTopology.prune_copies(repro_treefile, repro_binary_treefile)


                ro.r("library(phangorn)")
                ro.r("library(ape)")
                ro.r("expected_tree <- read.tree('{}')".format(expected_subsample_treefile))
                ro.r("reconstruct_tree <- read.tree('{}')".format(repro_binary_treefile))

                # FastTree makes unrooted tree; its root is arbitrarily set.  Reroot trees with same root before compare.
                # Use First tip in fasttree as common root.
                repro_tree = Phylo.read(repro_binary_treefile, 'newick')
                commonroot = repro_tree.get_terminals()[0].name

                ro.r("reconstruct_tree_rooted <- root(reconstruct_tree, outgroup='{}', resolve.root=TRUE)".format(commonroot))
                ro.r("expected_tree_rooted <- root(expected_tree, outgroup='{}', resolve.root=TRUE)".format(commonroot))
                #ro.r("png('reconstruct_tree_rooted.png')")
                #ro.r("plot(reconstruct_tree_rooted)")
                #ro.r("dev.off()")
                #ro.r("png('expected_tree_rooted.png')")
                #ro.r("plot(expected_tree_rooted)")
                #ro.r("dev.off()")

                #print ro.r("reconstruct_tree_rooted$tip.label[!reconstruct_tree_rooted$tip.label %in% expected_tree_rooted$tip.label]")
                #print ro.r("expected_tree_rooted$tip.label[!expected_tree_rooted$tip.label %in% reconstruct_tree_rooted$tip.label]")

                # robinson foulds distance
                ro.r("rf_dist <- RF.dist(expected_tree_rooted, reconstruct_tree_rooted)")
                rf_dist = ro.r("rf_dist")[0]



                self.assertEqual(rf_dist, 0,
                             "Reconstructed FastTree tree " + repro_binary_treefile + " should have same topology as input tree " +
                             expected_subsample_treefile + ". Instead robinson foulds distance =" + str(rf_dist))



    def test_subsample_final_popn(self):
        """
        Test that FastTree reconstruction of subsampled final population has same topology as
        the initial FastTree of the full final population (created with no topology constraint).
        """
        import rpy2.robjects as ro

        seed = self.configparser.getint(SECTION, "SEED")
        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.fasta"
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#

        final_popn_unconstrained_treefile = final_popn_tree_file.replace(".nwk", ".free.nwk")
        fasttree.fasttree_handler.make_tree(fasta_fname=final_popn_fasta,
                                            out_tree_fname=final_popn_unconstrained_treefile,
                                            threads=fasttree_threads, debug=True)

        for sample_fraction in [0.5]:
            subsample_fasta = final_popn_fasta.replace(".fasta", ".prune.{:.1}.fasta".format(sample_fraction))
            expected_subsample_treefile = final_popn_unconstrained_treefile.replace(".nwk", ".prune.{:.1}.nwk".format(sample_fraction))

            TestTopology.subsample(treefile=final_popn_unconstrained_treefile, fastafile=final_popn_fasta,
                                   out_treefile=expected_subsample_treefile, out_fastafile=subsample_fasta,
                                   seed=seed, fraction=sample_fraction)


            # Try to match the topology of the pruned final population topology when we use pruned fasta as input
            repro_treefile = fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                                 threads=fasttree_threads, debug=True)

            ro.r("library(phangorn)")
            ro.r("library(ape)")
            ro.r("expected_tree <- read.tree('{}')".format(expected_subsample_treefile))
            ro.r("reconstruct_tree <- read.tree('{}')".format(repro_treefile))

            # FastTree makes unrooted tree; its root is arbitrarily set.  Reroot trees with same root before compare.
            # Use first sequence in subsampled fasta as the root for both trees.
            first = SeqIO.parse(subsample_fasta, "fasta").next()

            ro.r("reconstruct_tree_rooted <- root(reconstruct_tree, outgroup='{}', resolve.root=TRUE)".format(first.id))
            ro.r("expected_tree_rooted <- root(expected_tree, outgroup='{}', resolve.root=TRUE)".format(first.id))

            # robinson foulds distance
            ro.r("rf_dist <- RF.dist(expected_tree_rooted, reconstruct_tree_rooted)")
            rf_dist = ro.r("rf_dist")[0]


            self.assertEqual(rf_dist, 0,
                         "Reconstructed FastTree tree " + repro_treefile + " should have same topology as input tree " +
                         expected_subsample_treefile + ". Instead robinson foulds distance =" + str(rf_dist))



    def test_resample_final_popn(self):
        """
        Test that FastTree reconstruction of subsampled with replacement final population has same topology as
        the initial FastTree of the pruned full final population (created with no topology constraint).
        """
        import rpy2.robjects as ro

        seed = self.configparser.getint(SECTION, "SEED")
        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.fasta"
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#

        final_popn_unconstrained_treefile = final_popn_tree_file.replace(".nwk", ".free.nwk")
        fasttree.fasttree_handler.make_tree(fasta_fname=final_popn_fasta,
                                            out_tree_fname=final_popn_unconstrained_treefile,
                                            threads=fasttree_threads, debug=True)

        for sample_fraction in [0.5]:
            subsample_fasta = final_popn_fasta.replace(".fasta", ".resample.{:.1}.fasta".format(sample_fraction))
            expected_subsample_treefile = final_popn_unconstrained_treefile.replace(".nwk", ".resample.{:.1}.nwk".format(sample_fraction))

            TestTopology.subsample(treefile=final_popn_unconstrained_treefile, fastafile=final_popn_fasta,
                                   out_treefile=expected_subsample_treefile, out_fastafile=subsample_fasta,
                                   seed=seed, fraction=sample_fraction)


            # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
            repro_treefile = fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                                 threads=fasttree_threads, debug=True)

            # robinson foulds distance only well defined for binary trees.  Prune sister copies of the same sequence from tree.
            # If the copies of the same sequence aren't polytomies in the tree, then FastTree screwed up.
            repro_binary_treefile =  repro_treefile.replace(".nwk", ".binary.nwk")
            TestTopology.prune_copies(repro_treefile, repro_binary_treefile)

            ro.r("library(phangorn)")
            ro.r("library(ape)")
            ro.r("expected_tree <- read.tree('{}')".format(expected_subsample_treefile))
            ro.r("reconstruct_tree <- read.tree('{}')".format(repro_binary_treefile))

            # FastTree makes unrooted tree; its root is arbitrarily set.  Reroot trees with same root before compare.
            # Use First tip in fasttree as common root.
            repro_tree = Phylo.read(repro_binary_treefile, 'newick')
            commonroot = repro_tree.get_terminals()[0].name

            ro.r("reconstruct_tree_rooted <- root(reconstruct_tree, outgroup='{}', resolve.root=TRUE)".format(commonroot))
            ro.r("expected_tree_rooted <- root(expected_tree, outgroup='{}', resolve.root=TRUE)".format(commonroot))
            #ro.r("png('reconstruct_tree_rooted.png')")
            #ro.r("plot(reconstruct_tree_rooted)")
            #ro.r("dev.off()")
            #ro.r("png('expected_tree_rooted.png')")
            #ro.r("plot(expected_tree_rooted)")
            #ro.r("dev.off()")

            #print ro.r("reconstruct_tree_rooted$tip.label[!reconstruct_tree_rooted$tip.label %in% expected_tree_rooted$tip.label]")
            #print ro.r("expected_tree_rooted$tip.label[!expected_tree_rooted$tip.label %in% reconstruct_tree_rooted$tip.label]")

            # robinson foulds distance
            ro.r("rf_dist <- RF.dist(expected_tree_rooted, reconstruct_tree_rooted)")
            rf_dist = ro.r("rf_dist")[0]



            self.assertEqual(rf_dist, 0,
                         "Reconstructed FastTree tree " + repro_binary_treefile + " should have same topology as input tree " +
                         expected_subsample_treefile + ". Instead robinson foulds distance =" + str(rf_dist))




if __name__ == "__main__":
    unittest.main()
