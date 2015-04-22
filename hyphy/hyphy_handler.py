import logging
import os
import subprocess
import config.settings as settings

LOGGER = logging.getLogger(__name__)

SELECTION_BF = "QuickSelectionDetection.bf"
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
    :param int threads:  HyPhy threads.  Note that HyPhy won't use more than 4 CPU for nucleotide model fitting and dN/dS calculations even if you give it more.
            Default 1.
    :param bool debug:  If True, then outputs hyphy stdout to log.
    :return:  path to HyPhy tab-separated output for sitewise dN/dS values.
        Each row in the tab-separate output represents a codon site in the alignment.
    :rtype: str
    """
    hyphy_filename_prefix = os.path.splitext(codon_fasta_filename)[0]  # Remove .fasta suffix
    hyphy_modelfit_filename = hyphy_filename_prefix + ".nucmodelfit"
    hyphy_dnds_tsv_filename = hyphy_filename_prefix + ".dnds.tsv"
    if debug:
        hyphy_log = hyphy_filename_prefix + ".hyphy.log"
    else:
        hyphy_log = os.devnull
    LOGGER.debug("Start HyPhy dN/dS " + hyphy_dnds_tsv_filename)
    if not os.path.exists(hyphy_dnds_tsv_filename) or os.path.getsize(hyphy_dnds_tsv_filename) <= 0:
        hyphy_input_str = "\n".join(["1",  # Universal
                                     "1",  # New analysis
                                     os.path.abspath(codon_fasta_filename),  # codon fasta
                                     "2",  #(2):[Custom] Use any reversible nucleotide model crossed with MG94.
                                     "012345",  # GTR
                                     os.path.abspath(tree_filename),  # tree file
                                     os.path.abspath(hyphy_modelfit_filename),  # model fit output file
                                     "3",  #(3):[Estimate] Estimate from data with branch corrections(slower).
                                     "1",  # Single Ancestor Counting
                                     "1",  # Full tree
                                     "1",  # Averaged
                                     "1",  # Approximate extended binomial distro
                                     "0.05",  # pvalue threshold for determining statistically significant dN/dS > 1 or dN/dS < 1.
                                     "2",  # Export to file
                                     os.path.abspath(hyphy_dnds_tsv_filename),  # dN/dS tsv output file
                                     "1\n"])  # Do not count approximate numbers of dN, dS rate classes supported by data

        # Feed window tree into hyphy to find dnds for the window
        with open(hyphy_log, 'w') as hyphy_log_fh:
            hyphy_cmd = [hyphy_exe, "BASEPATH=" + hyphy_basedir, "CPU=" + str(threads), SELECTION_BF]
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


def count_site_branch_subs(codon_fasta_filename, rooted_treefile, hyphy_exe=settings.DEFAULT_HYPHY_EXE, hyphy_basedir=settings.DEFAULT_HYPHY_BASEDIR,
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
    :param int threads:  HyPhy threads.  Note that HyPhy won't use more than 4 CPU for nucleotide model fitting and dN/dS calculations even if you give it more.
            Default 1.
    :param bool debug:  If True, then outputs hyphy stdout to log.
    :return:  path to HyPhy tab-separated output for sitewise dN/dS values.
        Each row in the tab-separate output represents a codon site in the alignment.
    :rtype: str
    """
    rooted_treefile_suffix = os.path.splitext(rooted_treefile)[1]
    nucmodelfit = rooted_treefile.replace(rooted_treefile_suffix, ".nucmodelfit")
    codonmodelfit = rooted_treefile.replace(rooted_treefile_suffix, ".codonmodelfit")
    leaf_anc_fasta = rooted_treefile.replace(rooted_treefile_suffix, ".anc.fasta")
    subst_tsv = rooted_treefile.replace(rooted_treefile_suffix, ".subst.tsv")

    # The hyphy alignment nexus file seems to always use the original tree instead of the whatever modifications by hyphy
    # Unfortunately, hyphy doesn't output branch lengths into the newick.  Instead, we have to get branch lengths from separate csvs
    # Check if branch lengths are actually constrained by hyphy's updated tree
    codon_treefile = leaf_anc_fasta.replace(".anc.fasta", ".codon.nwk")
    codontree_branchlen_csv = codon_treefile.replace(".nwk", ".branchlen.csv")

    LOGGER.debug("Counting site-branch subs for " + rooted_treefile)
    if (os.path.exists(leaf_anc_fasta) and os.path.getsize(leaf_anc_fasta) and
            os.path.exists(codon_treefile) and os.path.getsize(codon_treefile) and
            os.path.exists(codontree_branchlen_csv) and os.path.getsize(codontree_branchlen_csv)):
        LOGGER.warn("Not regenerating {} , {} , {}".format(leaf_anc_fasta,codon_treefile,  codontree_branchlen_csv))
    else:
        hyphy_input_str = "\n".join([codon_fasta_filename,  # msa aligned fasta of both samples
                                     rooted_treefile ,  # richard's root to tip rooted tree file.  length in nucleotide substitutions/site
                                     nucmodelfit,  # nucleotide model fit
                                     codonmodelfit, # codon model fit
                                     leaf_anc_fasta,  #  leaf and ancestor fasta
                                     subst_tsv,  # substitutions at each site-branch in tsv format
                                     codon_treefile,  # output tree in scaled to codon substitution/site lengths
                                     codontree_branchlen_csv  # csv for codon tree branch name and lengths

        ]) + "\n"

        # Feed window tree into hyphy to find ancestor sequences for the window
        hyphy_log = nucmodelfit.replace(".nucmodelfit", ".hyphy.log")
        with open(hyphy_log, 'w') as hyphy_log_fh:
            hyphy_cmd = [hyphy_exe, "CPU={}".format(threads),
                         "BASEPATH=" + hyphy_basedir,
                         "ReconstructAncestorSeq.bf"]
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


