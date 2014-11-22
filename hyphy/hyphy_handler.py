import logging
import sys
import os
import subprocess

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)

SELECTION_BF = "QuickSelectionDetection.bf"
NUC_MODEL_CMP_BS = "GTRrate.bf"
HYPHY_EXE = "HYPHYMP"
HYPHY_BASEDIR = "/usr/local/lib/hyphy/TemplateBatchFiles/"


# TODO:  multiple correction for pvalue
def calc_dnds(codon_fasta_filename, tree_filename, pvalue, hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEDIR, threads=1):
        hyphy_filename_prefix = os.path.splitext(codon_fasta_filename)[0]  # Remove .fasta suffix
        hyphy_modelfit_filename = hyphy_filename_prefix + ".nucmodelfit"
        hyphy_dnds_tsv_filename = hyphy_filename_prefix + ".dnds.tsv"
        hyphy_log = hyphy_filename_prefix + ".hyphy.log"
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
                                         "1",  # Single Acnestor Counting
                                         "1",  # Full tree
                                         "1",  # Averaged
                                         "1",  # Approximate extended binomial distro
                                         str(pvalue),  # pvalue threshold
                                         "2",  # Export to file
                                         os.path.abspath(hyphy_dnds_tsv_filename),  # dN/dS tsv output file
                                         "2\n"])  # Count approximate numbers of dN, dS rate classes supported by data

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


def calc_nuc_subst(codon_fasta_filename, tree_filename, hyphy_exe=HYPHY_EXE, hyphy_basedir=HYPHY_BASEDIR, threads=1):
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


