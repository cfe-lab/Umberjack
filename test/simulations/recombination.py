"""
Helper methods for handling recombination

If we need to make variations of the same tree to simulate recombination,
then for B breakpoints, make B+1 trees.
At each breakpoint, the population's genomes will follow a random different tree.

"""
import random
import config.settings
import logging
import math
import Bio.AlignIO
import Bio.Phylo


config.settings.setup_logging()
LOGGER = logging.getLogger(__name__)

MAX_SEED = math.pow(2, 32)-1  # internally asg_driver.py used numpy random generator which can only take up to 32bit seeds

def choose_breakpoints(genome_codons, num_breaks, seed=None):
    """
    Choose the locations of breakpoints within the genome.
    We define breakpoint as the 0-based codon position of the strand switch in recombination.
    First position in genome can't be breakpoint because there are no positions before it to switch strands from.
    Only allow breakpoints at codon start position (with respect to nucleotides) because
    we are simulating codon model in INDELible, which means that genome partitions can only begin at codon start.

    :param int genome_codons: number of codons in full genome
    :param int num_breaks:  number of recombination breakpoints to randomly select
    :param int seed:  random seed
    :return list:  list of integer tuples (start position, end position) of 1-based nucleotide positions of the
        recombination partitions.  Each partition represents a contiguous stretch on the same RNA strand.
        Recombination partitions will be listed in order.  Each new partition indicates that
        the opposite strand has been selected during recombination.
    """
    if seed is None:
        seed = random.randint(0, MAX_SEED)

    randomizer = random.Random(seed)
    LOGGER.info("Randomly selecting breakpoints with seed " + str(seed))


    # We define breakpoint as the 0-based codon start position of the strand switch in recombination.
    # First position in genome can't be breakpoint because there are no positions before it to switch strands from.
    breakpt_codons = sorted(randomizer.sample(range(1, genome_codons), num_breaks))


    # Go through the breakpoints and append the previous section into the list
    # Although breakpoints are given in 0-based codon positions,
    # we follow the Umbjerack custom of specifying genome partitions within filenames with 1-based nucleotide positions
    sections = []
    prev_sectn_start_nuc_base1 = 1
    for b in xrange(num_breaks):
        # Breakpoint indicates the start of another strand, aka another section
        breakpt_codon_base0 = breakpt_codons[b]
        sectn_start_nuc_base1 = (breakpt_codon_base0 * 3) + 1

        prev_sectn_end_nuc_base1 = sectn_start_nuc_base1 - 1
        sections.append((prev_sectn_start_nuc_base1, prev_sectn_end_nuc_base1))

        prev_sectn_start_nuc_base1 = sectn_start_nuc_base1

    prev_sectn_end_nuc_base1 = genome_codons * 3  # Append the section that starts from the last breakpoint and ends at genome-end
    sections.append((prev_sectn_start_nuc_base1, prev_sectn_end_nuc_base1))

    return sections


def get_sections_from_breakpoints(breakpoints, genome_codons):
    """
    Returns the section nuc 1based start and ends given the breakpoint 0-based codon positions.
    :param list breakpoints:  list of integer breakpoints
    :param int genome_codons:  total genome codons
    :return:
    """
    sections = []
    prev_sectn_start_nuc_base1 = 1
    for b in xrange(len(breakpoints)):
        # Breakpoint indicates the start of another strand, aka another section
        breakpt_codon_base0 = breakpoints[b]
        sectn_start_nuc_base1 = (breakpt_codon_base0 * 3) + 1

        prev_sectn_end_nuc_base1 = sectn_start_nuc_base1 - 1
        sections.append((prev_sectn_start_nuc_base1, prev_sectn_end_nuc_base1))

        prev_sectn_start_nuc_base1 = sectn_start_nuc_base1

    prev_sectn_end_nuc_base1 = genome_codons * 3  # Append the section that starts from the last breakpoint and ends at genome-end
    sections.append((prev_sectn_start_nuc_base1, prev_sectn_end_nuc_base1))

    return sections


def get_section_sizes(nuc_sections, is_codon_size=False):
    """
    Gets the sizes of each recombination section in nucleotide bases.
    :param list nuc_sections:  list of tuples (start nuc pos 1based, start end pos 1based) for each section
    :param bool is_codon_size:  whether to return the size in codons or nucleotides
    :return list:  integer list of sizes of each recombination section
    """
    sizes = []
    for start_nuc_base1, end_nuc_base1 in nuc_sections:
        nuc_size = end_nuc_base1 - start_nuc_base1 + 1
        if is_codon_size:
            codon_size = nuc_size / 3
            sizes.extend([codon_size])
        else:
            sizes.extend([nuc_size])
    return sizes



def make_recombo_tree(in_treefile, out_treefile=None, seed=None):
    """
    Simulates topology from recombination.
    Randomly selects a clade to prune and regraph to another randomly selected position in the tree.

    :param str in_treefile:  filepath to newick tree to be manipulated
    :param str out_treefile:  filepath to newick tree file to write to if given.  Uses Bio.Phylo to write tree to newick file.
    :param int seed:  random seed
    :return Bio.Phylo.TreeMixin:  pruned and regrafted clade, revised tree
    """
    if not seed:
        seed = random.randint(0, MAX_SEED)
    LOGGER.info("Randomizing prune regraph with seed " + str(seed))

    # Make sure that we can reproduce the tree by using a separate Random instance so we can seed it fresh,
    # without affecting other randomizations
    randomizer = random.Random(seed)

    # Randomly select a clade to move to another location in the tree.
    # Don't select the root clade (or the clade with all the tips), since it pretty much contains the entire tree.
    tree = Bio.Phylo.read(in_treefile, "newick")
    total_nodes = len(list(tree.find_clades()))
    relocate_clade_idx = randomizer.randint(1, total_nodes-1)  # in level (breadthfirst) order, the first node is the root node

    relocate_clade = prune_by_idx(tree=tree, idx=relocate_clade_idx, order="level")

    # Randomly select a remaining node to become the new sister of the pruned clade
    total_remaining_nodes = len(list(tree.find_clades()))
    # in level (breadthfirst) order, the first node is the root node
    # Don't regraft the pruned clade as root because that would mean
    # that the pruned clade is a recombinant of its previous parent and an ancestor even earlier than the previous root.
    dest_sister_idx = randomizer.randint(1, total_remaining_nodes-1)
    dest_sister_clade = list(tree.find_clades(order="level"))[dest_sister_idx]
    dest_new_par_br_len = randomizer.random() * dest_sister_clade.branch_length

    graft(tree=tree, src_subtree=relocate_clade, dest_sister_clade=dest_sister_clade, dest_new_par_br_len=dest_new_par_br_len)

    if out_treefile:
        Bio.Phylo.write(tree, out_treefile, "newick")

    return relocate_clade, tree


def prune_by_idx(tree, idx, order):
    """
    Prunes the clade with the 0-based index according to the tree-traversal order specified.
    Essentially copied and pasted from Bio.Phylo.BaseTree, but the BaseTree.prune() version only allows
    pruning of tips not clades.

    :param Bio.Phylo.BaseTree tree:  tree instance
    :param int idx:  0-based clade index
    :param str order:  tree traversal order.  [preorder, postorder, level]
    :return Bio.Phylo.Clade:  pruned clade
    """
    pruned_clade = None
    for i, clade in enumerate(tree.find_clades(order=order)):
        if i == idx:
            pruned_clade = clade
            path = tree.get_path(clade)  # path from (root, clade].  Ie exclude root, but include clade.
            if len(path) < 2:
                parent = path[-1]  # parent is root
            else:
                parent = path[-2]

            parent.clades.remove(clade)

            # Modified from Bio.Phylo.BaseTree.prune()
            # Adjust branch length of parent if pruning of clade makes parent a singleton node
            if len(parent.clades) == 1:
                # We deleted a branch from a bifurcation
                if parent == tree.root:
                    # If we're at the root, sister clade becomes new root
                    old_root_brlen = tree.root.branch_length
                    tree.root = parent.clades[0]
                    tree.root.branch_length += old_root_brlen
                else:
                    # If we're not at the root, collapse this parent
                    child = parent.clades[0]
                    if child.branch_length is not None:
                        child.branch_length += (parent.branch_length or 0.0)

                    if len(path) <= 2:
                        grandparent = tree.root
                    else:
                        grandparent = path[-3]

                    # Replace parent with child at the same place in grandparent
                    parent_index = grandparent.clades.index(parent)
                    grandparent.clades.pop(parent_index)
                    grandparent.clades.insert(parent_index, child)

            break

    return pruned_clade


def graft(tree, src_subtree, dest_sister_clade, dest_new_par_br_len):
    """
    Grafts a subree as a sister to the destination clade.
    If the destination sister clade is an inner node, inserts an intermediate node in between destination sister clade
    and its parent.
    New intermediate node will have branch length of dest_br_len.  Shortens the destination sister clade branch length
    so that the total tree length is the same.
    If the new intermediate node zero branch, simply graphs the subtree
    as a polytomy sister of the destination sister clade without creating an intermediate node.

            src_clade:
           |===
        ===|
           |===

        tree:  I = location of new intermediate node

            |--I--tip1
        --N0|
            |---tip2

        Graft src_clade as a sister clade to tip1.  Insert new intermediate node in between tip1 and its parent N0
        Result:

               |--tip1
               |
            |--I   |===
            |  |===|
            |      |===
        --N0|
            |---tip2



    :param Bio.Phylo.TreeMixin tree:  tree instance to be grafted onto
    :param Bio.Phylo.TreeMixin src_subtree:  subtree to graft onto tree
    :param Bio.Phylo.Clade dest_sister_clade:  clade that will become sister to grafted clade
    :param float dest_new_par_br_len:  branch length of new intermediate node.  Must be within [0, branch length of dest_sister_clade)
    :return Bio.Phylo:  revised tree instance
    :raises ValueError: if the new intermediate branch length >= destination sister clade branch length,
        or if the destination sister clade is the root node.
    """

    if dest_sister_clade == tree.root:
        raise ValueError("Grafting to root branch not allowed")
    if dest_new_par_br_len >= dest_sister_clade.branch_length:
        raise ValueError("New intermediate node not allowed to be on top or past destination sister node. "
                         "Instead specify the child of dest_sister_clade and " +
                         "specify dest_br_len in between [0, branch length of that child)")
    if dest_new_par_br_len < 0:
        raise ValueError("New intermediate node branch length must be within [0, dest_sister_clade branch length)")

    # Case:  new intermediate node branch length = 0
    # Result:  Graft src_clade as polytomy sister to dest_sister_clade, don't bother creating intermediate node


    # Case:  new intermediate node branch length  in (0, dest_sister_clade branch length)
    # Result:  Insert new intermediate node between dest_sister_clade and its parent.
    #   Graft src_clade as child of new intermediate node and sister to dest_sister_clade.

    # path between (root--> dest_sister_clade], excluding root including dest_sister_clade
    node_path = tree.get_path(dest_sister_clade)
    if len(node_path) >= 2:
        parent = node_path[-2]
    else:
        parent = tree.root  #  parent is root


    # new intermediate node is directly on top of the parent of the destination sister clade
    if dest_new_par_br_len == 0:
        # Don't bother inserting a new node in between the parent and the current clade
        # because it would just have zero branch length.
        # Make subtree a polytomy sister of destination sister clade
        parent.clades.append(src_subtree)
    # new intermediate node is in between parent and destination sister clade
    else:
        # Create a new node in between destination sister clade and its parent.
        # New node will become parent of subtree and destination sister clade
        parent.split(n=1, branch_length=dest_new_par_br_len)
        newnode = parent.clades[-1]
        parent.clades.remove(dest_sister_clade)
        newnode.clades = [src_subtree.root, dest_sister_clade]
        # Keep total branch length from destination sister clade to its previous parent the same
        dest_sister_clade.branch_length -= dest_new_par_br_len

    return tree


def write_section_fasta(full_popn_fasta, nuc_sections, filename_prefix):
    """
    Create fastas for each recombinant section
    :param str full_popn_fasta: filepath to full population fasta for every individual in population
    :param list nuc_sections:  list of (start nuc 1-based position, end nuc 1based position) of contiguous sections of DNA
        separated by recombinant breakpoints
    :param str filename_prefix:  filepath and filename prefix of output fasta for each recombinant section.
        Will be appended with <nuc start>_<nuc_end>.fasta for each recombinant section
    :return list:  list of filepaths to output fastas
    """
    output_fastas = []
    aln =  Bio.AlignIO.read(full_popn_fasta, "fasta")
    for section_start_base1, section_end_base1 in nuc_sections:  # 1-based nucleotide recombinant section positions
        # Slicing Bio.AlignIO alignments uses 0-based nucleotide coordinates
        full_popn_section_fasta = (filename_prefix + ".{}_{}.fasta".format(section_start_base1, section_end_base1))
        Bio.AlignIO.write(aln[:, section_start_base1-1:section_end_base1], full_popn_section_fasta, "fasta")
        output_fastas.append(full_popn_section_fasta)

    return output_fastas
