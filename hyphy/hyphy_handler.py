import logging
import os
import subprocess
import config.settings as settings

LOGGER = logging.getLogger(__name__)

SELECTION_BF = "QuickSelectionDetection.bf"
COUNT_SUBS_BF = "CountSubs.bf"
NUC_MODEL_CMP_BS = "GTRrate.bf"


# Columns in the HyPhy dN/dS tab-separates values file
HYPHY_TSV_DN_COL = 'dN'
HYPHY_TSV_DS_COL = 'dS'
HYPHY_TSV_S_COL = 'Observed S Changes'
HYPHY_TSV_N_COL = 'Observed NS Changes'
HYPHY_TSV_SCALED_DN_MINUS_DS_COL = 'Scaled dN-dS'
HYPHY_TSV_PROB_FULLSEQ_S_COL = 'P{S leq. observed}'
HYPHY_TSV_NEG_PROB_FULLSEQ_S_COL = 'P{S geq. observed}'
HYPHY_TSV_EXP_S_COL = 'E[S Sites]'
HYPHY_TSV_EXP_N_COL = 'E[NS Sites]'


def calc_dnds(codon_fasta_filename, tree_filename, hyphy_exe=settings.DEFAULT_HYPHY_EXE, hyphy_basedir=settings.DEFAULT_HYPHY_BASEDIR,
              hyphy_libdir=settings.DEFAULT_HYPHY_LIBDIR,
              threads=1, debug=False):
    """
    Calculates sitewise dN/dS using Single-most Likely Ancestor Counting with branch corrections.
    Averages substitutions across all possible codons for ambiguous codons.
    Fits a general time reversible nucleotide model crossed with MG94 codon model and saves it in the same directory as the fasta file with suffix ".nucmodelfit".
    Redirects HyPhy stdout/stderr to the same directory as the codon fasta file with the suffix ".hyphy.log".
    HyPhy outputs sitewise dN/dS in tab-separated file in the same directory as the fasta file with the suffix ".dnds.tsv"
    Note that if there are stop codons at any site, HyPhy will quietly remove those sites from the output.
    It is best to mask the stop codons with "NNN" before feeding them into this method.

    Will not overwrite existing HyPhy output for the same fasta file.

    :param str codon_fasta_filename:  path to multiple sequence aligned fasta file that starts at the beginning of a codon.
    :param str tree_filename: path to phylogenetic tree for sequences in fasta.
    :param str hyphy_exe:  path to multithreaded (but not MPI enabled) hyphy executable HYPHYMP.  If empty, uses HYPHYMP from PATH.
    :param str hyphy_basedir:  path to HyPhy TemplateBatchFiles directory.  If empty, uses /usr/local/lib/hyphy/TemplateBatchFiles
    :param str hyphy_libdir:  path to HyPhy Custom LIBDIR containing the out of box hyphy batch files.
    :param int threads:  HyPhy threads.  Note that HyPhy won't use more than 4 CPU for nucleotide model fitting and dN/dS calculations even if you give it more.
            Default 1.
    :param bool debug:  If True, then outputs hyphy stdout to log.
    :return:  path to HyPhy tab-separated output for sitewise dN/dS values.
        Each row in the tab-separate output represents a codon site in the alignment.
    :rtype: str
    """
    hyphy_filename_prefix = os.path.splitext(codon_fasta_filename)[0]  # Remove .fasta suffix
    hyphy_modelfit_filename = hyphy_filename_prefix + ".codonmodelfit"
    hyphy_dnds_tsv_filename = hyphy_filename_prefix + ".dnds.tsv"
    hyphy_codon_br_len_tsv_filename = hyphy_filename_prefix + ".codon.brlen.csv"
    hyphy_anc_fasta = hyphy_filename_prefix + ".anc.fasta"

    if debug:
        hyphy_log = hyphy_filename_prefix + ".hyphy.log"
    else:
        hyphy_log = os.devnull
    LOGGER.debug("Start HyPhy dN/dS " + hyphy_dnds_tsv_filename)
    if not os.path.exists(hyphy_dnds_tsv_filename) or os.path.getsize(hyphy_dnds_tsv_filename) <= 0:
        hyphy_input_str = "\n".join([ "-1",  # No per-site-branch substitution TSV file
                                     os.path.abspath(hyphy_dnds_tsv_filename),  # dN/dS tsv output file
                                    "1",  # Universal Genetic Code
                                     os.path.abspath(codon_fasta_filename),  # codon fasta
                                     os.path.abspath(tree_filename),  # tree file
                                     os.path.abspath(hyphy_modelfit_filename),  # model fit output file
                                     os.path.abspath(hyphy_codon_br_len_tsv_filename), # codon tree branch length csv
                                     os.path.abspath(hyphy_anc_fasta),  # reconstructed ancestor + tips fasta
                                     "\n"])

        # Feed window tree into hyphy to find dnds for the window
        with open(hyphy_log, 'w') as hyphy_log_fh:
            hyphy_cmd = [hyphy_exe,
                         "BASEPATH=" + hyphy_basedir,
                         "LIBPATH=" + hyphy_libdir,
                         "CPU=" + str(threads),
                         COUNT_SUBS_BF]
            hyphy_proc = subprocess.Popen(hyphy_cmd, stdin=subprocess.PIPE, stdout=hyphy_log_fh,
                                          stderr=hyphy_log_fh,
                                          shell=False, env=os.environ)
            hyphy_proc.communicate(hyphy_input_str)

            if hyphy_proc.returncode:
                raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)

        LOGGER.debug("Done HyPhy dN/dS " + hyphy_dnds_tsv_filename)
    else:
        LOGGER.debug("Found existing HyPhy for window " + hyphy_dnds_tsv_filename + ". Not regenerating")

    return hyphy_dnds_tsv_filename


def count_site_branch_subs(codon_fasta_filename, rooted_treefile,
                           hyphy_exe=settings.DEFAULT_HYPHY_EXE, hyphy_basedir=settings.DEFAULT_HYPHY_BASEDIR,
                           hyphy_libdir=settings.DEFAULT_HYPHY_LIBDIR,
                           threads=1, debug=False):
    """
    Calculates sitewise dN/dS using Single-most Likely Ancestor Counting with branch corrections.
    Averages substitutions across all possible codons for ambiguous codons.
    Fits a general time reversible nucleotide model crossed with MG94 codon model and saves it in the same directory as the fasta file with suffix ".nucmodelfit".
    Redirects HyPhy stdout/stderr to the same directory as the codon fasta file with the suffix ".hyphy.log".
    HyPhy outputs sitewise dN/dS in tab-separated file in the same directory as the fasta file with the suffix ".dnds.tsv"
    Note that if there are stop codons at any site, HyPhy will quietly remove those sites from the output.
    It is best to mask the stop codons with "NNN" before feeding them into this method.

    Will not overwrite existing HyPhy output for the same fasta file.

    :param str codon_fasta_filename:  path to multiple sequence aligned fasta file that starts at the beginning of a codon.
    :param str rooted_treefile: path to phylogenetic tree for sequences in fasta.
    :param str hyphy_exe:  path to multithreaded (but not MPI enabled) hyphy executable HYPHYMP.  If empty, uses HYPHYMP from PATH.
    :param str hyphy_basedir:  path to HyPhy Custom Batchfile Dir.
    :param str hyphy_libdir:  path to HyPhy Custom LIBDIR containing the out of box hyphy batch files.
    :param int threads:  HyPhy threads.  Note that HyPhy won't use more than 4 CPU for nucleotide model fitting and dN/dS calculations even if you give it more.
            Default 1.
    :param bool debug:  If True, then outputs hyphy stdout to log.
    :return:  path to HyPhy tab-separated output for sitewise dN/dS values.
        Each row in the tab-separate output represents a codon site in the alignment.
    :rtype: str
    """
    rooted_treefile_suffix = os.path.splitext(rooted_treefile)[1]
    codonmodelfit = rooted_treefile.replace(rooted_treefile_suffix, ".codonmodelfit")
    leaf_anc_fasta = rooted_treefile.replace(rooted_treefile_suffix, ".anc.fasta")
    subst_tsv = rooted_treefile.replace(rooted_treefile_suffix, ".subst.tsv")

    # The hyphy alignment nexus file seems to always use the original tree instead of the whatever modifications by hyphy
    # Unfortunately, hyphy doesn't output branch lengths into the newick.  Instead, we have to get branch lengths from separate csvs
    # Check if branch lengths are actually constrained by hyphy's updated tree
    codontree_branchlen_csv = rooted_treefile_suffix.replace(rooted_treefile_suffix, ".codon.brlen.csv")

    LOGGER.debug("Counting site-branch subs for " + rooted_treefile)
    if (os.path.exists(leaf_anc_fasta) and os.path.getsize(leaf_anc_fasta) and
            os.path.exists(codontree_branchlen_csv) and os.path.getsize(codontree_branchlen_csv) and
            os.path.exists(subst_tsv) and os.path.getsize(subst_tsv)):
        LOGGER.warn("Not regenerating {} , {} , {}".format(leaf_anc_fasta, codontree_branchlen_csv, subst_tsv))
    else:
        if debug:
            hyphy_log = rooted_treefile.replace(rooted_treefile_suffix, ".hyphy.log")
        else:
            hyphy_log = os.devnull

        LOGGER.debug("Start HyPhy Site-Branch Sub Count " + subst_tsv)

        hyphy_input_str = "\n".join([ subst_tsv,  # Per-site-branch substitution TSV file
                                      "-1",  # No dN/dS tsv output file
                                      "1",  # Universal Genetic Code
                                      os.path.abspath(codon_fasta_filename),  # codon fasta
                                      os.path.abspath(rooted_treefile),  # tree file
                                      os.path.abspath(codonmodelfit),  # model fit output file
                                      os.path.abspath(codontree_branchlen_csv), # codon tree branch length csv
                                      os.path.abspath(leaf_anc_fasta),  # reconstructed ancestor + tips fasta
                                      "\n"])

        # Feed window tree into hyphy to find dnds for the window
        with open(hyphy_log, 'w') as hyphy_log_fh:
            hyphy_cmd = [hyphy_exe,
                         "BASEPATH=" + hyphy_basedir,
                         "LIBPATH=" + hyphy_libdir,
                         "CPU=" + str(threads),
                         COUNT_SUBS_BF]
            hyphy_proc = subprocess.Popen(hyphy_cmd, stdin=subprocess.PIPE, stdout=hyphy_log_fh,
                                          stderr=hyphy_log_fh,
                                          shell=False, env=os.environ)
            hyphy_proc.communicate(hyphy_input_str)

            if hyphy_proc.returncode:
                raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)

        LOGGER.debug("Done Counting site-branch subs for  " + rooted_treefile)

    return subst_tsv


def calc_nuc_subst(codon_fasta_filename, tree_filename, hyphy_exe=settings.DEFAULT_HYPHY_EXE, hyphy_basedir=settings.DEFAULT_HYPHY_BASEDIR, threads=1):
    """
    Fits a general time reversible (4 discrete gamma rates) nucleotide model and saves it in same directory as fasta with suffix ".nucmodelfit".
    :param str codon_fasta_filename:  path to multiple sequence aligned fasta file that starts at the beginning of a codon.
    :param str tree_filename: path to phylogenetic tree for sequences in fasta.
    :param str hyphy_exe:  path to multithreaded (but not MPI enabled) hyphy executable HYPHYMP.  If empty, uses HYPHYMP from PATH.
    :param str hyphy_basedir:  path to HyPhy TemplateBatchFiles directory.  If empty, uses /usr/local/lib/hyphy/TemplateBatchFiles
    :param int threads:  HyPhy threads.  Note that HyPhy won't use more than 4 CPU for nucleotide model fitting and dN/dS calculations even if you give it more.
            Default 1.
    """
    hyphy_filename_prefix = os.path.splitext(codon_fasta_filename)[0]  # Remove .fasta suffix
    hyphy_modelfit_filename = hyphy_filename_prefix + ".nucmodelfit"
    hyphy_log = hyphy_filename_prefix + ".hyphy.log"
    LOGGER.debug("Start HyPhy nucleotide model " + hyphy_modelfit_filename)
    if not os.path.exists(hyphy_modelfit_filename) or os.path.getsize(hyphy_modelfit_filename) <= 0:
        hyphy_input_str = "\n".join([
            os.path.abspath(codon_fasta_filename),  # codon fasta
            os.path.abspath(tree_filename),  # tree file
            "4",  # Number of rate classes in rate variation models (e.g. 4):  # TODO:  make this configurable
            os.path.abspath(hyphy_modelfit_filename),  # model fit output file
            "\n"])

        with open(hyphy_log, 'w') as hyphy_log_fh:
            hyphy_cmd = [hyphy_exe, "BASEPATH=" + hyphy_basedir, "CPU=" + str(threads), NUC_MODEL_CMP_BS]
            hyphy_proc = subprocess.Popen(hyphy_cmd, stdin=subprocess.PIPE, stdout=hyphy_log_fh,
                                          stderr=hyphy_log_fh,
                                          shell=False, env=os.environ)
            hyphy_proc.communicate(hyphy_input_str)

            if hyphy_proc.returncode:
                raise subprocess.CalledProcessError(cmd=hyphy_cmd, returncode=hyphy_proc.returncode)
        LOGGER.debug("Done HyPhy nucleotide model " + hyphy_modelfit_filename)
    else:
        LOGGER.debug("Found existing HyPhy nucleotide model " + hyphy_modelfit_filename + ". Not regenerating")


