
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
import logging
import Utility
import math
import dendropy  # calculates weighted robinson foulds distance
import dendropy.calculate.treecompare

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

LOGGER = logging.getLogger(__name__)

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

    @staticmethod
    def get_rf_dist(expected_treefile, actual_treefile, is_reroot=False):
        """

        :param expected_treefile:
        :param actual_treefile:
        :param is_reroot:
        :return:
        """
        # Use R phangorn package to calculate
        import rpy2.robjects as ro
        ro.r("library(phangorn)")
        ro.r("library(ape)")
        ro.r("expected_tree <- read.tree('{}')".format(expected_treefile))
        ro.r("actual_tree <- read.tree('{}')".format(actual_treefile))

        # Check if trees have different labels.  R will throw error if they do.
        ro.r("expected_notin_actual <- expected_tree$tip.label[!expected_tree$tip.label %in% actual_tree$tip.label]")
        expected_notin_actual = ro.r("expected_notin_actual")
        rlen_e_notin_a = ro.r("length(expected_notin_actual)")[0]
        if rlen_e_notin_a !=  0:
            raise ValueError(
                             "There are tips in expected tree not in actual tree=" + str(expected_notin_actual) +
                             ". expected tree = " + expected_treefile +
                             " actual tree=" + actual_treefile)

        ro.r("actual_notin_expected <- actual_tree$tip.label[!actual_tree$tip.label %in% expected_tree$tip.label]")
        actual_notin_expected = ro.r("actual_notin_expected")
        rlen_a_notin_e = ro.r("length(actual_notin_expected)")[0]
        if rlen_a_notin_e != 0:
            raise ValueError(
                             "There are tips in actual tree not in expected tree=" + str(actual_notin_expected) +
                             ".  expected tree = " + expected_treefile +
                             " actual tree=" + actual_treefile)

        if is_reroot:
            ro.r("newroot <- expected_tree$tip.label[1]")
            ro.r("expected_tree <- root(expected_tree, outgroup=newroot, resolve.root=TRUE)")
            ro.r("actual_tree <- root(actual_tree, outgroup=newroot, resolve.root=TRUE)")


        # robinson foulds distance
        ro.r("rf_dist <- RF.dist(expected_tree, actual_tree)")
        #print ro.r("rf_dist")
        rf_dist = ro.r("rf_dist")[0]
        return rf_dist


    @staticmethod
    def get_weighted_rf_dist(expected_treefile, actual_treefile):
        """
        Used dendropy to find robinson foulds distance weighted by the branch length difference between trees
        of the edge that defines a bifurcation in a tree.
        If the bifurcation does not exist in the other tree, then the other tree's branch length is set to zero.
        Allows for polytomies and unroots the trees by default before calculating bifurcations.

        :param expected_treefile:
        :param actual_treefile:
        :return:
        """
        tns = dendropy.TaxonNamespace()  # both trees need to have the same taxons

        expected_tree = dendropy.Tree.get(
                path=expected_treefile,
                schema='newick',
                taxon_namespace=tns)
        actual_tree = dendropy.Tree.get(
                path=actual_treefile,
                schema='newick',
                taxon_namespace=tns)

        rf_dist = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(expected_tree, actual_tree)

        return rf_dist

    @staticmethod
    def get_weighted_rf_dist_from_str(expected_treestr, actual_treestr):
        """
        Used dendropy to find robinson foulds distance weighted by the branch length difference between trees
        of the edge that defines a bifurcation in a tree.
        If the bifurcation does not exist in the other tree, then the other tree's branch length is set to zero.
        Allows for polytomies and unroots the trees by default before calculating bifurcations.

        :param str expected_treestr: newick string
        :param str actual_treestr: newick string
        :return:
        """
        tns = dendropy.TaxonNamespace()  # both trees need to have the same taxons

        expected_tree = dendropy.Tree.get(
                data=expected_treestr,
                schema='newick',
                taxon_namespace=tns)
        actual_tree = dendropy.Tree.get(
                data=actual_treestr,
                schema='newick',
                taxon_namespace=tns)

        rf_dist = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(expected_tree, actual_tree)

        return rf_dist


    def cmp_topology(self, expected_treefile, actual_treefile, is_reroot=False):
        """
        Compares the topology of the trees.  Asserts if there is non-zero Robinson Foulds distance.
        If reroot specifed, then uses the first tip in the expected tree as the root.
        :return:
        """
        rf_dist = TestTopology.get_rf_dist(expected_treefile, actual_treefile, is_reroot)

        self.assertEqual(rf_dist, 0,
                     "Actual tree " + actual_treefile + " should have same topology as expected tree " +
                     expected_treefile + ". Instead robinson foulds distance =" + str(rf_dist))





    @staticmethod
    def intersect_size(range_start1, range_end1, range_start2, range_end2):
        """
        :param int range_start1:  first range start position
        :param int range_end1:  first range end position
        :param int range_start2:  2nd range start position
        :param int range_end2:  2nd range end position
        :return int : Size of the intersection or negative number of they don't intersect
        """
        return min(range_end1, range_end2) - max(range_start1, range_start2)


    @staticmethod
    def intersect_range(range_start1, range_end1, range_start2, range_end2):
        """
        :param int range_start1:  first range start position
        :param int range_end1:  first range end position
        :param int range_start2:  2nd range start position
        :param int range_end2:  2nd range end position
        :return int, int : start and end position of the intersection
        """
        return max(range_start1, range_start2), min(range_end1, range_end2)


    @staticmethod
    def resample_fasta_from_seq(src_fasta, dest_fasta, out_fasta, slice_range=None):
        """
        Takes the sequence names from src_fasta.
        Then resamples the sequences in dest_fasta so that they correspond to src_fasta.
        Slices dest_fasta according to the slice coordinates.
        Writes the resampled dest_fasta to out_fasta.

        Assumes that
        - the names in src_fasta follow format [base name]_read[num].
        - the names in dest_fasta follow format [base name]

        :param str src_fasta:   fasta containing the sequence names to copy from
        :param str dest_fasta: fasta containing the sequence names to modify to match src_fasta
        :param str out_fasta:   output fasta to write modified dest_fasta
        :param tuple slice_range:  (int 1-based start wrt dest_fasta, int 1-based end wrt fasta)
        If not supplied, then entire dest_fasta written to out_fasta.
        """
        with open(out_fasta, 'w') as fh_out:
            dest_recs = SeqIO.to_dict(SeqIO.parse(dest_fasta, "fasta"))
            for src_rec in SeqIO.parse(src_fasta, "fasta"):
                rec_basename = src_rec.id.split("_")[0]
                fh_out.write(">" + src_rec.id + "\n")

                dest_rec = dest_recs.get(rec_basename)

                if slice_range:
                    slice_start_base1, slice_end_base1 = slice_range
                    fh_out.write(str(dest_rec.seq[slice_start_base1-1:slice_end_base1]) + "\n")
                else:
                    fh_out.write(str(dest_rec.seq) + "\n")


    @staticmethod
    def calc_window_tree_dist(sim_data, window_fasta, window_treefile,
                              win_start, win_end):
        """
        Calculates the distance between a window tree and the recombinant trees in the sections that overlap the window
        as the weighted average robinson folds - branchlength distance.
        :param SimData sim_data:  simulated data instance
        :param str window_fasta: filepath to umberjack window fasta
        :param str window_treefile: filepath to umberjack window tree
        :param int win_start:  1 based nucleotide start position of window
        :param int win_end:  1 based nucleotide end position of window
        :return float:  weighted average robinson foulds-branchlength distance between the window tree
        and all full population recombinant trees that overlap the window
        """

        full_popn_fasta = sim_data.get_fasta()

        # Make a temp directory in the umberjack output directory to hold the expected trees for the recombinant section-window overlap
        expected_dir = os.path.dirname(window_treefile) + os.sep + "expected"
        if not os.path.exists(expected_dir):
            os.makedirs(expected_dir)

        expected_prefix = expected_dir + os.sep + os.path.splitext(os.path.basename(full_popn_fasta))[0]

        weighted_ave_dist = 0.0
        for full_popn_recombo_treefile in sim_data.get_recombo_trees():
            recomb_start, recomb_end = sim_data.get_recombo_breaks_from_file(full_popn_recombo_treefile)

            if recomb_start > win_end or recomb_end < win_start:
                continue

            recomb_win_start, recomb_win_end = TestTopology.intersect_range(recomb_start, recomb_end, win_start, win_end)

            # Full population Treefile with tips resampled according to the read sequences used in the window fasta
            # Branch lengths can't be trusted if the recombinant-window overlap sequences
            # are different from the recombinant sequences.  In which case, only used as for topology guide.
            expected_topology_treefile = expected_prefix + ".expected.topology.{}_{}.nwk".format(recomb_win_start, recomb_win_end)
            # Full population Treefile with tips resampled according to the read sequences used in the window fasta
            # Branch lengths rescaled to fit the sequences in recombinant section-window overlap
            expected_treefile = expected_prefix + ".expected.{}_{}.nwk".format(recomb_win_start, recomb_win_end)

            # For each read in the window, make a copy of the individual's tip in the full population tree
            TestTopology.resample_tree_from_fasta(treefile=full_popn_recombo_treefile,
                                                  out_treefile=expected_topology_treefile,
                                                  fastafile=window_fasta)

            # If the recombinant-window overlap is the same as the recombinant section,
            # then we can trust the recombinant-window branch lengths
            recombo_win_size = recomb_win_end - recomb_win_start + 1
            win_size = win_end - win_start + 1
            if recomb_win_start == recomb_start and recomb_win_end == recomb_end:
                shutil.move(expected_topology_treefile, expected_treefile)
            else:

                # For each read in the window, make a copy of the individual's sequence from the full population fasta
                expected_fasta = expected_prefix + ".expected.{}_{}.fasta".format(recomb_win_start, recomb_win_end)
                TestTopology.resample_fasta_from_seq(src_fasta=window_fasta, dest_fasta=full_popn_fasta,
                                                     out_fasta=expected_fasta,
                                                     slice_range=(recomb_win_start, recomb_win_end))

                # Recalculate branch lengths but constrain topology
                # since the recombinant-window overlap sequences different from the recombinant section sequences
                fasttree.fasttree_handler.make_tree_repro(fasta_fname=expected_fasta,
                                             intree_fname=expected_topology_treefile,
                                             outtree_fname=expected_treefile)

            rf_dist = TestTopology.get_weighted_rf_dist(expected_treefile, window_treefile)

            weighted_ave_dist +=  (recombo_win_size/float(win_size)) * rf_dist

        # Clean up if we haven't encountered any errors.  Leave the expected files there for debugging if we crash before this.
        if os.path.exists(expected_dir):
            shutil.rmtree(expected_dir)


        return weighted_ave_dist





    def test_full_popn_tree(self):
        """
        Tests that the FastTree tree of the full population uses the same topology as the input tree generated from the coalescent simulator.
        :return:
        """

        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.fasta"
        # /data/umberjack_unittest/mixed/umberjack_unittest.mixed.nwk
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#
        nodup_final_popn_tree_file = final_popn_tree_file.replace(".nwk", ".nodup.nwk")

        # /data/umberjack_unittest/umberjack_unittest.rename.nwk
        orig_tree_file = SIM_DATA_DIR + os.sep + SIM_DATA_FILENAME_PREFIX + ".rename.nwk"
        nodup_orig_tree_file = orig_tree_file.replace(".nwk", ".nodup.nwk")

        # INDELible will output identical sequences when branch lengths are too short, even if the original coalescent tree
        # specifies non-zero branch lengths.  As a result, FastTree puts these identical sequences as polytomies.
        # Get rid of tips from duplicated INDELible sequences from both the original coalescent and final population tree.
        TestTopology.prune_copies_by_seq(fastafile=final_popn_fasta, in_treefile=final_popn_tree_file, out_treefile=nodup_final_popn_tree_file)
        TestTopology.prune_copies_by_seq(fastafile=final_popn_fasta, in_treefile=orig_tree_file, out_treefile=nodup_orig_tree_file)

        self.cmp_topology(expected_treefile=nodup_orig_tree_file, actual_treefile=nodup_final_popn_tree_file)



    def test_full_popn_free_tree(self):
        """
        Tests that the FastTree tree of the full population made without topology constraints
        uses the same topology as the input tree generated from the coalescent simulator.
        :return:
        """
        # /data/umberjack_unittest/mixed/umberjack_unittest.mixed.nwk
        final_popn_free_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.free.nwk"
        nodup_final_popn_free_tree_file = final_popn_free_tree_file.replace(".nwk", ".nodup.nwk")
        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.fasta"

        # /data/umberjack_unittest/umberjack_unittest.rename.nwk
        orig_tree_file = SIM_DATA_DIR + os.sep + SIM_DATA_FILENAME_PREFIX + ".rename.nwk"
        nodup_orig_tree_file = orig_tree_file.replace(".nwk", ".nodup.nwk")

        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        fasttree.fasttree_handler.make_tree(fasta_fname=final_popn_fasta,
                                                             out_tree_fname=final_popn_free_tree_file,
                                                             threads=fasttree_threads, debug=True)

        # INDELible will output identical sequences when branch lengths are too short, even if the original coalescent tree
        # specifies non-zero branch lengths.  As a result, FastTree puts these identical sequences as polytomies.
        # Get rid of tips from duplicated INDELible sequences from both the original coalescent and final population tree.
        TestTopology.prune_copies_by_seq(fastafile=final_popn_fasta, in_treefile=final_popn_free_tree_file, out_treefile=nodup_final_popn_free_tree_file)
        TestTopology.prune_copies_by_seq(fastafile=final_popn_fasta, in_treefile=orig_tree_file, out_treefile=nodup_orig_tree_file)

        self.cmp_topology(expected_treefile=nodup_orig_tree_file, actual_treefile=nodup_final_popn_free_tree_file)



    def test_indelible_topology(self):
        """
        Test that the INDELible scaled tree has the same topology as the original coalescent tree
        :return:
        """
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

            self.cmp_topology(expected_treefile=orig_tree_file, actual_treefile=indelible_nwk)



    def test_fasttree_topology(self):
        """
        Test that FastTree can reproduce the same topology of an INDELible tree given the INDELible fasta sequences as input
        and the INDELIble tree as a constraint.
        We do this on the full genome population scaled at only 1 mutation scaling rate
        so that we don't have multiple mutation rates per sequence confounding issues with tree reconstruction.
        :return:
        """

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

            # At low mutation rates, INDELible will output fasta sequences that are identical,
            # and won't fully correspond to the input tree topology.
            # Remove the duplicate copies before comparing topologies.
            nodup_indelible_nwk = indelible_nwk.replace(".nwk", ".nodup.nwk")
            nodup_repro_treefile = repro_treefile.replace(".nwk", ".nodup.nwk")
            TestTopology.prune_copies_by_seq(fastafile=indelible_tip_fasta, in_treefile=indelible_nwk, out_treefile=nodup_indelible_nwk)
            TestTopology.prune_copies_by_seq(fastafile=indelible_tip_fasta, in_treefile=repro_treefile, out_treefile=nodup_repro_treefile)

            self.cmp_topology(expected_treefile=nodup_indelible_nwk, actual_treefile=nodup_repro_treefile)


    @staticmethod
    def miss_seq(seq, fraction, seed=None, is_rand_fraction=False):
        """
        Gaps-out or N's out fraction of seq.
        Randomly assigns gaps at end, N's in the middle.
        :param str seq:  nucleotide sequence
        :param float fraction:  fraction of bases to turn into gaps or Ns.  If is_rand_fraction, then randomly assigns UP TO fraction of bases.
        Otherwise, assigns floor (fraction of bases).
        :param bool is_rand_fraction:  whether to randomize the fraction of bases to assign as gaps or N's.
        :return str:  gapped, N'ed sequence.
        """
        randomizer = random.Random(seed)
        if is_rand_fraction:
            max_miss = math.floor(fraction * len(seq))
            total_miss = randomizer.randint(0, max_miss)
        else:
            total_miss = int(math.floor(fraction * len(seq)))

        total_miss_end = randomizer.randint(0, total_miss)
        total_left_miss = randomizer.randint(0, total_miss_end)
        total_right_miss =  total_miss_end - total_left_miss

        total_mid_miss = total_miss - total_miss_end
        mid_size = len(seq) - total_miss_end

        # print total_left_miss
        # print total_right_miss
        # print total_mid_miss

        miss_seq = "-" * total_left_miss

        mid_miss_idx = randomizer.sample(range(total_left_miss, total_left_miss+mid_size), total_mid_miss)
        for i in range(total_left_miss, total_left_miss+mid_size):
            if i in mid_miss_idx:
                miss_seq += "N"
            else:
                miss_seq += seq[i]

        miss_seq += "-" * total_right_miss

        return miss_seq


    @staticmethod
    def subsample(treefile, fastafile, out_treefile, out_fastafile, sample_fraction, replace=False, seed=None,
                  miss_fraction=None, removedup=True):
        """
        Subsamples population with or without replacement up to the desired fraction.  Returns subsampled tree and fasta.
        Removes duplicated sequences from the output tree (ie sequences in the original fasta that are duplicated, or sequences that have been resampled).
        Does not remove duplicated sequences from the output fasta.
        :param str treefile:  filepath to newick population tree
        :param str fastafile:  filepath to population sequences in fasta format
        :param float sample_fraction:  fraction of population to sample
        :param bool replace:  whether to sample with replacement
        :param int seed:  the random seed to randomly select individuals.
        :param float miss_fraction:  fraction of sequence to turn into gaps or N's.
        :param bool removedup:  whether to remove duplicated sequences from the output tree.
            NB:  Duplicated sequences will always be output in output fasta.
        :return [str]:  list of sequence names selected. Renames sequences to  "<name>_read<copy">
        """
        # NB:  ART read simulator generates read with names like otu1-read2.  Use underscore instead of hyphen to make HyPhy friendly.
        tree = Phylo.read(treefile, "newick")
        tips = tree.get_terminals()

        sample_size = int(round(sample_fraction *  len(tips)))

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
                for i in range(0, count):
                    new_tipname = "{}_read{}".format(tipname, i)
                    tipname_sample.extend([new_tipname])
                    fh_out_fasta.write(">" + new_tipname + "\n")

                    seq = str(seqdict[tipname].seq)
                    if miss_fraction:  # Gap-out, N-out bases to simulate missing data
                        miss_seq = TestTopology.miss_seq(seq=seq, fraction=miss_fraction)
                    else:
                        miss_seq = seq
                    fh_out_fasta.write(miss_seq + "\n")



        # Prune the tree of tips that weren't selected for the sample.
        for tip in tips:
            if tip.name not in selecttip_to_count.keys() or selecttip_to_count[tip.name] == 0:
                tree.prune(tip)

        # Append resampled tips to the tree as polytomies.
        # Rename tips to [tipname]_read[index]
        for tip in tree.get_terminals():
            dup_count = selecttip_to_count[tip.name]
            if dup_count >= 2:
                tip.split(n=dup_count, branch_length=0)
                dup_tips = tip.clades
            else:
                dup_tips = [tip]
            for i, dup_tip in enumerate(dup_tips):
                dup_tip.name = "{}_read{}".format(tip.name, i)

        if removedup:
            # We expected the resampled tips to be identical to each other, resulting in polytomies in the output tree.
            # But if we replace the tip fasta sequences with random missing data, then copied sequences
            # will no longer be identical.
            # Check which sequences are actually identical, then remove the copies from the output tree,
            # keeping the first copy encountered (in alphabetical tip name order).
            dup_ids = TestTopology.find_dup_seq(fasta=out_fastafile)
            for tip in tree.get_terminals():
                if tip.name in dup_ids:
                    tree.prune(tip)


        Phylo.write(tree, out_treefile, "newick")

        return tipname_sample

    @staticmethod
    def resample_tree_from_fasta(treefile, fastafile, out_treefile):
        """
        Fasta contains resampled tips, renamed with format [tipname]_read[num].
        Resample the tips in the tree according to the resampled seq in the fasta.
        Write out to out_treefile
        :param str treefile:  filepath to newick population tree
        :param str fastafile:  filepath to population sequences in fasta format
        :param str out_treefile:  filepath to output newick tree
        """
        # NB:  ART read simulator generates read with names like otu1-read2.  Use underscore instead of hyphen to make HyPhy friendly.
        tree = Phylo.read(treefile, "newick")
        tips = tree.get_terminals()

        tipbasename_to_newtips = {}
        for rec in SeqIO.parse(fastafile, "fasta"):
            tip_basename = rec.id.split("_")[0]
            if tip_basename not in tipbasename_to_newtips:
                tipbasename_to_newtips[tip_basename] = []
            tipbasename_to_newtips[tip_basename].extend([rec.id])

        # Prune the tree of tips that weren't selected for the sample.
        for tip in tips:
            if tip.name not in tipbasename_to_newtips:
                tree.prune(tip)

        # Append resampled tips to the tree as polytomies.
        # Rename tips to [tipname]_read[index]
        for tip in tree.get_terminals():
            if tip.name not in tipbasename_to_newtips:
                raise ValueError("tip should be in tree")
            dup_count = len(tipbasename_to_newtips[tip.name])
            if dup_count >= 2:
                tip.split(n=dup_count, branch_length=0)
                dup_tips = tip.clades
            else:
                dup_tips = [tip]
            for i, dup_tip in enumerate(dup_tips):
                newtip_name = tipbasename_to_newtips[tip.name][i]
                dup_tip.name = newtip_name

            # If we split new tips and the previous tip has zero branch length, then collapse it
            if dup_count >= 2 and tip.branch_length == 0:
                tree.collapse(tip)

        Phylo.write(tree, out_treefile, "newick")




    def test_subsample(self):
        """
        Test that FastTree reconstruction of subsampled INDELible population has same topology as INDELible tree.
        Don't use final population whose sequences are concatenated from multiple mutation rates.
        Instead, just use a population whose entire genome has been scaled by a single mutation rate.

        """

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
                                       sample_fraction=sample_fraction, seed=seed)



                # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
                repro_treefile = fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                                     threads=fasttree_threads, debug=True)
                # Get rid of duplicate sequences in the tips that are a result of low mutation rate.
                # These show up as polytomies in the actual tree, even if the expected tree has non-zero branch lengths between them.
                nodup_repro_treefile = repro_treefile.replace(".nwk", ".nodup.nwk")
                TestTopology.prune_copies_by_seq(in_treefile=repro_treefile, out_treefile=nodup_repro_treefile, fastafile=subsample_fasta)
                self.cmp_topology(expected_treefile=expected_subsample_treefile, actual_treefile=nodup_repro_treefile, is_reroot=False)



    @staticmethod
    def find_dup_seq(fasta):
        """
        Finds duplicated sequences (headers can be different, but sequences are the same) in a fasta.
        Uses the first copy (ordered by alphabetical order by header id) as the original copy.
        :param fasta:
        :return [str]:  header ids of sequences that are duplicates of original seq.
        """
        rec_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        seq_set = set([])
        dup_ids = []
        for id in sorted(rec_dict.keys()):
            seq = str(rec_dict[id].seq)
            if seq in seq_set:
                dup_ids.extend([id])
            else:
                seq_set.add(seq)

        return dup_ids



    @staticmethod
    def prune_copies_by_tipname(treefile, out_treefile):
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


    @staticmethod
    def prune_copies_by_seq(in_treefile, out_treefile, fastafile):
        """
        This is a helper function to account for the fact that if there is insufficient branch length,
        INDELible will output multiple identical sequences, even if the desired input tree has non-zero branches lengths.

        If there are identical sequences in the actual fasta file, then they will result in polytomies in the actual tree.
        We need to ensure that the expected_treefile is revised to remove the extra copies of the identical sequences.
        Only the first copy of the sequence in alpha order by tip name is kept.

        :return :
        """
        dup_ids = TestTopology.find_dup_seq(fasta=fastafile)
        if dup_ids:
            LOGGER.warn("Found duplicate sequences in fasta " + fastafile)

        tree = Phylo.read(in_treefile, 'newick')
        for tip in tree.get_terminals():
            if tip.name in dup_ids:
                tree.prune(tip)

        Phylo.write(tree, out_treefile, "newick")



    def test_subsample_replace(self):
        """
        Test that FastTree reconstruction of subsampled with replacement INDELible population has same topology as INDELible tree
        that has been pruned for unsampled tips.
        Don't use final population whose sequences are concatenated from multiple muitation rates.
        Instead, just use a population whose entire genome has been scaled by a single mutation rate.
        :return:
        """
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
                                       sample_fraction=sample_fraction, replace=True, seed=seed)



                # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
                repro_treefile = fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                                     threads=fasttree_threads, debug=True)

                # Get rid of duplicate sequences in the tips that are a result of low mutation rate.
                # These show up as polytomies in the actual tree, even if the expected tree has non-zero branch lengths between them.
                nodup_repro_treefile =  repro_treefile.replace(".nwk", ".nodup.nwk")
                TestTopology.prune_copies_by_seq(in_treefile=repro_treefile, out_treefile=nodup_repro_treefile, fastafile=subsample_fasta)


                self.cmp_topology(expected_treefile=expected_subsample_treefile, actual_treefile=nodup_repro_treefile, is_reroot=False)



    def test_subsample_final_popn(self):
        """
        Test that FastTree reconstruction of subsampled final population has same topology as
        the initial FastTree of the full final population (created with no topology constraint).
        """
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
                                   sample_fraction=sample_fraction, seed=seed, replace=False, removedup=False)


            # Try to match the topology of the pruned final population topology when we use pruned fasta as input
            repro_treefile = subsample_fasta.replace(".fasta", ".repro.nwk")
            fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                out_tree_fname=repro_treefile,
                                                threads=fasttree_threads, debug=True)

            self.cmp_topology(expected_treefile=expected_subsample_treefile, actual_treefile=repro_treefile, is_reroot=False)



    def test_resample_final_popn(self):
        """
        Test that FastTree reconstruction of subsampled with replacement final population has same topology as
        the initial FastTree of the pruned full final population (created with no topology constraint).
        """
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
                                   sample_fraction=sample_fraction, seed=seed, replace=True, removedup=False)


            # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
            repro_treefile = subsample_fasta.replace(".fasta", ".repro.nwk")
            fasttree.fasttree_handler.make_tree(fasta_fname=subsample_fasta,
                                                out_tree_fname=repro_treefile,
                                                threads=fasttree_threads, debug=True)
            self.cmp_topology(expected_treefile=expected_subsample_treefile, actual_treefile=repro_treefile, is_reroot=False)



    def test_missing_seq_final_popn(self):
        """
        Checks that Ns or gaps do not drastically alter topology.
        :return:
        """
        seed = self.configparser.getint(SECTION, "SEED")
        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.fasta"
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#

        final_popn_unconstrained_treefile = final_popn_tree_file.replace(".nwk", ".free.nwk")
        fasttree.fasttree_handler.make_tree(fasta_fname=final_popn_fasta,
                                            out_tree_fname=final_popn_unconstrained_treefile,
                                            threads=fasttree_threads, debug=True)

        for missing_fraction in [0.125]:
            missing_fasta = final_popn_fasta.replace(".fasta", ".missing.{:.3}.fasta".format(missing_fraction))
            expected_subsample_treefile = final_popn_unconstrained_treefile.replace(".nwk", ".missing.{:.3}.nwk".format(missing_fraction))

            TestTopology.subsample(treefile=final_popn_unconstrained_treefile, fastafile=final_popn_fasta,
                                   out_treefile=expected_subsample_treefile, out_fastafile=missing_fasta,
                                   sample_fraction=1, seed=seed, miss_fraction=missing_fraction,
                                   replace=False, removedup=False)


            # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
            repro_treefile = missing_fasta.replace(".fasta", ".repro.nwk")
            fasttree.fasttree_handler.make_tree(fasta_fname=missing_fasta,
                                                out_tree_fname=repro_treefile,
                                                threads=fasttree_threads, debug=True)

            self.cmp_topology(expected_treefile=expected_subsample_treefile, actual_treefile=repro_treefile, is_reroot=False)


    def test_missing_seq_resample_final_popn(self):
        """
        Checks that random Ns or gaps in resampled sequences do not drastically alter topology.
        :return:
        """
        seed = self.configparser.getint(SECTION, "SEED")
        fasttree_threads = self.configparser.getint(SECTION, "PROCS")

        final_popn_fasta = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX + ".mixed.fasta"
        final_popn_tree_file = SIM_DATA_DIR + os.sep + "mixed" + os.sep + SIM_DATA_FILENAME_PREFIX  + ".mixed.nwk"#

        final_popn_unconstrained_treefile = final_popn_tree_file.replace(".nwk", ".free.nwk")
        fasttree.fasttree_handler.make_tree(fasta_fname=final_popn_fasta,
                                            out_tree_fname=final_popn_unconstrained_treefile,
                                            threads=fasttree_threads, debug=True)

        for sample_fraction in [0.5]:
            for missing_fraction in [0.125]:
                sub_fasta = final_popn_fasta.replace(".fasta", ".resample.{:.1}.missing.{:.3}.fasta".format(sample_fraction, missing_fraction))
                expected_subsample_treefile = final_popn_unconstrained_treefile.replace(".nwk",
                                                                                        ".resample.{:.1}.missing.{:.3}.nwk".format(sample_fraction, missing_fraction))

                TestTopology.subsample(treefile=final_popn_unconstrained_treefile, fastafile=final_popn_fasta,
                                       out_treefile=expected_subsample_treefile, out_fastafile=sub_fasta,
                                       sample_fraction=sample_fraction, seed=seed, miss_fraction=missing_fraction,
                                       replace=True, removedup=False)


                # Try to reproduce the INDELible tree with the INDELible fasta sequences as input
                repro_treefile = sub_fasta.replace(".fasta", ".repro.nwk")
                fasttree.fasttree_handler.make_tree(fasta_fname=sub_fasta,
                                                    out_tree_fname=repro_treefile,
                                                    threads=fasttree_threads, debug=True)
                self.cmp_topology(expected_treefile=expected_subsample_treefile, actual_treefile=repro_treefile, is_reroot=False)


if __name__ == "__main__":

    # for i in range(0, 10):
    #     print 'test ' + str(i)
    #     #print TestTopology.miss_seq("AAAAACCCCCTTTTTGGGGG", fraction=0.125, seed=i)
    #     print TestTopology.miss_seq("AAAAAAAAAACCCCCCCCCCTTTTTTTTTTGGGGGGGGGG", fraction=0.125, is_rand_fraction=True)

    unittest.main()
