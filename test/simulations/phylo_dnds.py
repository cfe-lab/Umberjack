# Calculates the true dn/ds along a phylogeny where all the leaf sequences, ancestral sequences, and phylogeny are known.
from Bio import Phylo
from Bio import SeqIO as SeqIO
from Bio import Seq as Seq
import itertools
import Utility



def calc_total_poss_subst(codon):
    """
    Calculate the total of possible point mutations that yield synonymous amino acid substitution.
    :param Bio.Seq.Seq codon:  3bp nucleotide codon
    :return int: total of point mutations that yield synonymous amino acid substitutions
    """
    orig_aa = Seq.translate(codon)
    total_poss_syn = 0
    total_poss_nonsyn = 0
    for codon_pos in range(0, Utility.NUC_PER_CODON):
        nuc = codon[codon_pos]
        for mut_str in ("A", "C", "T", "G"):
            mut = Seq.Seq(mut_str)
            if str(mut).upper() == str(nuc).upper():
                continue
            mut_codon = codon[:codon_pos] + mut + codon[codon_pos+1:]
            mut_aa = Seq.translate(mut_codon)
            if str(orig_aa).upper() == str(mut_aa).upper():
                total_poss_syn += 1
            else:
                total_poss_nonsyn += 1

    return total_poss_syn, total_poss_nonsyn


def calc_total_subst(codon1, codon2):
    """
    Returns total synonymous substitutions, nonsynonymous substitutions
    :param Bio.Seq.Seq codon1:  3bp codon
    :param Bio.Seq.Seq codon2:  3bp codon
    :return tuple (int, int):  (total point mutations that yield same amino acid, total point mutations that yield different amino acid)
    """
    total_syn = 0
    total_nonsyn = 0
    orig_aa1 = Seq.translate(codon1)
    for pos, nucstr1 in enumerate(str(codon1)):
        nucstr2 = str(codon2[pos])
        if nucstr1.upper() != nucstr2.upper():
            mut_codon1 = codon1[:pos] + Seq.Seq(nucstr2) + codon1[pos+1:]
            mut_aa1 = Seq.translate(mut_codon1)
            if str(orig_aa1) == str(mut_aa1):
                total_syn += 1
            else:
                total_nonsyn += 1
    return total_syn, total_nonsyn





def is_subst(codon1, codon2):
    """
    Returns if there is a substitution between pair of nucleotide codons
    :param Bio.Seq.Seq codon1:  3bp codon
    :param Bio.Seq.Seq codon2: 3bp codon
    :return bool: True if there is a substitution between the nucleotide codons.  False otherwise
    """
    return str(codon1) != str(codon2)



def get_parent(tree, child_clade):
    """
    Returns the parent clade of the given clade
    :param Bio.Phylo.Tree tree:
    :param Bio.Phylo.Clade child_clade:
    :return:
    """
    node_path = tree.get_path(child_clade)  # list of clades in path from just after root to the child clade
    return node_path[-2]  # the last item in list is the child clade.  The 2nd last item is the parent clade.


def calc_popn_dnds(leaf_fasta, ancestor_fasta, treefile):
    """
    Calculates the number of synonymous and nonsynonymous substitutions at each amino acid site for the population.
    :param str leaf_fasta: full filepath of fasta containing the  sequences of each individual in population (leaf nodes in tree)
    :param str ancestor_fasta: full filepath of fasta containing the sequences of each ancestor in the population (inner nodes in tree)
    :param file treefile: full filepath of newick tree file showing phylogeny of individuals and ancestors
    :return list of floats [float, float, float...]:  list of dn/ds for each amino acid site
    """
    tree = Phylo.read(treefile, "newick")
    nuc_recs = SeqIO.to_dict(itertools.chain(SeqIO.parse(leaf_fasta, "fasta"), SeqIO.parse(ancestor_fasta, "fasta")))
    num_codons = len(nuc_recs[nuc_recs.keys()[0]].seq)/Utility.NUC_PER_CODON

    # num_codons = 20  # TODO:  delete me
    all_site_dnds = [0.0]*num_codons
    all_site_total_poss_syn = [0.0]*num_codons
    all_site_total_poss_nonsyn = [0.0]*num_codons
    all_site_syn = [0.0]*num_codons
    all_site_nonsyn = [0.0]*num_codons

    # Traverse each branch (parent-child pair).
    # Count number of nonsynonymous substitutions, synonymous substitutions between each parent-child pair at each codon site.
    #root_seq = nuc_recs[tree.clade.name].seq
    #for leaf in tree.get_terminals():
        # path_clades = tree.get_path(leaf)  # list of clades in path from just below root to leaf
        # parent_seq = root_seq

    innernodes = tree.find_clades(terminal=False)
    for parent in innernodes:
        #print "Processing path=" + ",".join([clade.name for clade in path_clades])
        #for child in path_clades:  # traverse from all child clades from just below root to leaf

        parent_seq = nuc_recs[parent.name].seq
        for child in parent:
            print "Processing parent-child=" + parent.name + "-" + child.name
            child_seq = nuc_recs[child.name].seq

            for codon_site in range(0, num_codons):  #0-based codon sites
                nuc_pos = codon_site*Utility.NUC_PER_CODON
                parent_codon_seq = parent_seq[nuc_pos:nuc_pos + Utility.NUC_PER_CODON]
                child_codon_seq = child_seq[nuc_pos:nuc_pos + Utility.NUC_PER_CODON]

                total_codon_poss_syn, total_codon_poss_nonsyn = calc_total_poss_subst(parent_codon_seq)
                if is_subst(parent_codon_seq, child_codon_seq):
                    total_codon_syn, total_codon_nonsyn = calc_total_subst(parent_codon_seq, child_codon_seq)
                else:
                    total_codon_syn, total_codon_nonsyn = 0,0

                all_site_total_poss_syn[codon_site] += total_codon_poss_syn
                all_site_total_poss_nonsyn[codon_site] += total_codon_poss_nonsyn
                all_site_syn[codon_site] += total_codon_syn
                all_site_nonsyn[codon_site] += total_codon_nonsyn

                print "site=" + str(codon_site) + " par={} child={} poss_syn={} poss_nonsyn={} syn={} nonsyn={}".format(parent_codon_seq, child_codon_seq,
                                                                                                                        total_codon_poss_syn, total_codon_poss_nonsyn,
                                                                                                        total_codon_syn, total_codon_nonsyn)



    for codon_site in range(0, num_codons):
        if all_site_total_poss_nonsyn[codon_site] != 0 and all_site_total_poss_syn[codon_site] != 0:
            ave_site_nonsyn = all_site_nonsyn[codon_site]/all_site_total_poss_nonsyn[codon_site]
            ave_site_syn = all_site_syn[codon_site]/all_site_total_poss_syn[codon_site]
            if ave_site_syn != 0:
                all_site_dnds[codon_site] = ave_site_nonsyn/ave_site_syn

    return all_site_dnds
