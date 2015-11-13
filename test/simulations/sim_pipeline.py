"""
The full pipeline for generating simulated population reads for unit testing.
Usage:  python sim_pipeline.py [config file]
"""

import subprocess
import os
import logging
import sys
import ConfigParser
import hyphy.hyphy_handler as hyphy_handler
import fasttree.fasttree_handler as fasttree_handler
import config.settings
import indelible.indelible_handler
from collections import OrderedDict
import Utility
import recombination
import Bio.AlignIO

LOGGER = logging.getLogger(__name__)
SECTION = "sim"


def get_path_str(path, pardir):
    """
    If absolute path, then returns the path as is.
    If relative path, then returns absolute path of concatenated pardir/path
    :param str path:  absolute or relative file or directory path
    :param str pardir: parent directory to concatenate to path if path is relative directory
    :return str: absolute resolved path
    """
    if not os.path.isabs(path):
        return os.path.join(pardir, path)
    else:
        return path


def cwd():
    """
    :return str:  current directory
    """
    return os.path.dirname(os.path.realpath(__file__))


def call_cmd(cmd, prereqs=None, outputs=None):
    """
    Simulates makefile functionality.
    Checks if outputs don't exist or prereqs have been modified after outputs.
    If so, then executes the command to regenerate the outputs.
    :param [str] cmd:  commandline arguments to execute
    :param list prereqs:  list of prerequisite filepaths.  If any don't exist or have been updated later than any of the outputs,
        then executes the command to update the outputs.
    :param list outputs:  list of output filepaths.
    :return:
    """
    # Set PYTHONPATH for Umberjack test simulation scripts
    umberjack_dir = os.path.abspath(cwd() + os.sep + os.pardir + os.sep + os.pardir)
    test_sim_dir = os.path.abspath(umberjack_dir + os.sep + "test" + os.sep + "simulations")
    indelible_dir = os.path.abspath(umberjack_dir + os.sep + "test" + os.sep + "simulations" + os.sep + "indelible")
    os.environ["PYTHONPATH"] = ":".join([
        umberjack_dir,
        test_sim_dir,
        indelible_dir,
        "$PYTHONPATH"])

    latest_prereq_date = None
    if prereqs:
        for prereq in prereqs:
            if latest_prereq_date is None or os.path.getmtime(prereq) > latest_prereq_date:
                latest_prereq_date  = os.path.getmtime(prereq)

    is_up_to_date = True
    if outputs:
        for output in outputs:
            if not os.path.exists(output) or (latest_prereq_date is not None and latest_prereq_date > os.path.getmtime(output)):
                is_up_to_date = False
                break

    if is_up_to_date:
        LOGGER.warn("Not regenerating " + ",".join(outputs))
    else:
        LOGGER.debug("About to execute " + " ".join(cmd))
        subprocess.check_call(cmd, env=os.environ)
        LOGGER.debug("Finished execute " + " ".join(cmd))



if __name__ == "__main__":
    config.settings.setup_logging()  # Set up logging formatting using default Umberjack logging settings

    config_file = sys.argv[1]
    cfgparser = ConfigParser.RawConfigParser()
    cfgparser.read(config_file)

    OUTDIR = os.path.dirname(config_file)  # Output directory for simulated data
    SEED = cfgparser.getint(SECTION, "SEED")
    FILENAME_PREFIX = cfgparser.get(SECTION, "FILENAME_PREFIX")
    NUM_CODON_SITES = cfgparser.getint(SECTION, "NUM_CODON_SITES")
    NUM_INDIV = cfgparser.getint(SECTION, "NUM_INDIV")
    if cfgparser.has_option(SECTION, "NUM_BREAKPOINTS"):
        NUM_BREAKPOINTS = cfgparser.getint(SECTION, "NUM_BREAKPOINTS")
    else:
        NUM_BREAKPOINTS = 0


    OUT_TOPOLOGY_DIR = OUTDIR + os.sep + "topology"
    OUT_FULL_POPN_DIR =OUTDIR + os.sep + "fullpopn"
    OUT_READ_DIR = OUTDIR + os.sep + "reads"
    OUT_ALN_DIR = OUTDIR + os.sep + "aln"
    OUT_SUBS_DIR = OUTDIR + os.sep + "subs"


    #######################
    # Decide where recombination breakpoints should be
    ########################

    # 1-based nucleotide positions of each contiguous DNA strands separated by recombination breakpoint
    recombo_sections = recombination.choose_breakpoints(genome_codons=NUM_CODON_SITES, num_breaks=NUM_BREAKPOINTS, seed=SEED)


    #######################
    # Generate Tree Topology for Full Population
    ######################

    if not os.path.exists(OUT_TOPOLOGY_DIR):
        os.makedirs(OUT_TOPOLOGY_DIR)

    # For each recombinant section, create different phylogeny from ancestral selection graphs.
    # These form the input topology for INDELible simulated population sequences.
    # Note that INDELible only uses the input topology as a guide and may stray away from it.
    # An example is if the branch lengths are very short between tips, INDELible may just output identical tip sequences,
    # resulting in polytomies, whereas the input trees will always be binary trees.
    #
    # Relabel tree nodes to more manageable names.  Reformat tree so that indelible can handle it.
    # Filenames for the breakpoint trees contain the breakpoints, which are randomly selected.
    # We need to ask for the breakpoint tree filenames given a random seed.
    # We can't ask them to print it out to stdout because the LOGGER might log to stdout
    relabel_phylogeny_exe = os.path.abspath(cwd() + os.sep + "relabel_phylogeny.py")


    renamed_treefiles = []
    break_treefile_prefixes = []
    if not recombo_sections:
        break_treefile_prefixes = [OUT_TOPOLOGY_DIR + os.sep + FILENAME_PREFIX ]
        renamed_treefiles = [OUT_TOPOLOGY_DIR + os.sep + FILENAME_PREFIX  + ".rename.nwk"]

    for section_nuc_start_base1, section_nuc_end_base1 in recombo_sections:
        break_treefile_prefix = OUT_TOPOLOGY_DIR + os.sep + FILENAME_PREFIX + ".break.{}_{}".format(section_nuc_start_base1, section_nuc_end_base1)
        break_treefile_prefixes.extend([break_treefile_prefix])
        renamed_treefiles.extend([break_treefile_prefix + ".rename.nwk"])

    # Generate the breakpoint treefiles for the entire population if they are out of date
    relabel_phylogeny_cmd = ["python", relabel_phylogeny_exe,
                             "-f", ",".join(break_treefile_prefixes),
                             "-t", str(NUM_INDIV),
                             "-s", str(SEED)]
    call_cmd(cmd=relabel_phylogeny_cmd, prereqs=[], outputs=renamed_treefiles)


    #######################
    # Generate Genome Sequence for Full Population
    ######################

    # Use Indelible to create population sequences at different scaling factors (ie mutation rates)
    INDELIBLE_BIN_DIR = get_path_str(cfgparser.get(SECTION, "INDELIBLE_BIN_DIR"), OUTDIR)
    INDELIBLE_SCALING_RATES = cfgparser.get(SECTION, "INDELIBLE_SCALING_RATES")

    # For each recombinant section, divvy up again into different substitution rates.
    # Each substitution rate is assigned to a contiguous block of codons whose size
    # is drawn from a negative binomial with parameters  p=probability of success=1/(total distinct substitution rates),
    # and n=number of trials = total genome codons.
    # Each tree topology - substitution rate combination forms a distinct INDELible partition.



    indelible_partition_csv = OUT_FULL_POPN_DIR + os.sep +"partition.csv"
    # true alignment of entire extent population sequences
    full_popn_fasta = OUT_FULL_POPN_DIR + os.sep +"{}_TRUE.fasta".format(FILENAME_PREFIX)
    # MSA alignment of ancestors of extent population
    ancestor_fasta = OUT_FULL_POPN_DIR + os.sep +"{}_ANCESTRAL.fasta".format(FILENAME_PREFIX)
    # CSV with columns Site, Class, Partition, Is_Inserted, Omega, ScalingRate, TreeFile
    # describing the dn/ds, tree topology, mutation scaling rate used at each site
    rates_csv =  OUT_FULL_POPN_DIR + os.sep +"{}_RATES.csv".format(FILENAME_PREFIX)
    # Population consensus
    full_popn_cons_fasta = OUT_FULL_POPN_DIR + os.sep + FILENAME_PREFIX + ".consensus.fasta"

    if not os.path.exists(OUT_FULL_POPN_DIR):
        os.makedirs(OUT_FULL_POPN_DIR)

    if os.path.exists(indelible_partition_csv) and os.path.getsize(indelible_partition_csv):
        LOGGER.warn("Not regenerating " + indelible_partition_csv)
    else:
        partition_sizes = recombination.get_section_sizes(nuc_sections=recombo_sections, is_codon_size=True)
        treefile_to_codons = OrderedDict()
        for i, treefile in enumerate(renamed_treefiles):
            codons = partition_sizes[i]
            treefile_to_codons[treefile] = codons

        scaling_rates = [float(x) for x in INDELIBLE_SCALING_RATES.split(",")]
        indelible.indelible_handler.write_partition_csv(partition_csv=indelible_partition_csv,
                                                        treefile_to_codons=treefile_to_codons,
                                                        tree_scaling_rates=scaling_rates,
                                                        seed=SEED)

    batch_indelible_exe = os.path.abspath(cwd() + "/indelible/batch_indelible.py")
    indelible_cmd = ["python", batch_indelible_exe,
                     "-p", indelible_partition_csv, # Describes the genome partitions with different tree topology or mutation rates
                     "-s", str(SEED),  # random seed
                     "-f", FILENAME_PREFIX,  # Output Filename prefix without the directory
                     "-o", OUT_FULL_POPN_DIR,  # indelible output file directory
                     "-i", INDELIBLE_BIN_DIR]  # indelible bin dir

    call_cmd(cmd=indelible_cmd, prereqs=renamed_treefiles +  [indelible_partition_csv],
             outputs=[full_popn_fasta, ancestor_fasta, rates_csv, full_popn_cons_fasta])


    ########################################
    # Simulate MiSeq reads from the population genomes.
    # Align them against the full population consensus
    ###############################################

    generate_reads_exe = os.path.abspath(cwd() + os.sep + "generate_reads.py")
    generate_reads_cmd = ["python", generate_reads_exe,
                          config_file]
    reads_fq1 = OUT_READ_DIR + os.sep + FILENAME_PREFIX + ".reads.1.fq"
    reads_fq2 = OUT_READ_DIR+ os.sep + FILENAME_PREFIX + ".reads.2.fq"
    art_sam = OUT_READ_DIR + os.sep + FILENAME_PREFIX + ".reads.sam"
    art_errfree_sam = OUT_READ_DIR + os.sep + FILENAME_PREFIX + ".reads.errFree.sam"

    aln_sam = OUT_ALN_DIR + os.sep + FILENAME_PREFIX + ".reads.consensus.bwa.sort.query.sam"

    call_cmd(cmd=generate_reads_cmd, prereqs=[full_popn_fasta, ancestor_fasta, rates_csv],
             outputs=[reads_fq1, reads_fq2, art_sam, art_errfree_sam, full_popn_cons_fasta, aln_sam])


    #################################
    # Infer topology and site dN/dS from full population sequences actually generated from INDELible
    #################################

    # TODO:
    # Don't use FastTree to infer the true tree.
    # We have the ancestral sequences.  Just count the substitutions between each parent-child.  Divide by sites.
    # This is true branch length.
    if not os.path.exists(OUT_SUBS_DIR):
        os.makedirs(OUT_SUBS_DIR)


    FASTTREE_EXE = get_path_str(cfgparser.get(SECTION, "FASTTREE_EXE"), OUTDIR)

    PROCS = cfgparser.getint(SECTION, "PROCS")
    HYPHY_EXE = get_path_str(cfgparser.get(SECTION, "HYPHY_EXE"), OUTDIR)
    HYPHY_BASEPATH = get_path_str(cfgparser.get(SECTION, "HYPHY_BASEPATH"), OUTDIR)


    # FastTree requires fasta sequences for each tree
    full_popn_section_prefix = OUT_FULL_POPN_DIR + os.sep + FILENAME_PREFIX + ".break"
    full_popn_section_fastas = recombination.write_section_fasta(full_popn_fasta=full_popn_fasta,
                                                                nuc_sections=recombo_sections,
                                                                filename_prefix=full_popn_section_prefix)

    for i, (section_nuc_start_base1, section_nuc_end_base1) in enumerate(recombo_sections):
        full_popn_section_fasta = full_popn_section_fastas[i]
        renamed_treefile = renamed_treefiles[i]

        # Use FastTree to recalculate the branch lengths constrained by the same topology as the asg_driver.py topologies
        inferred_full_popn_section_treefile = (OUT_SUBS_DIR + os.sep +
                                       FILENAME_PREFIX + ".break.{}_{}.fasttree.nwk".format(section_nuc_start_base1, section_nuc_end_base1))

        sample_genomes_tree_fname = fasttree_handler.make_tree_repro(fasta_fname=full_popn_fasta,
                                                                     intree_fname=renamed_treefile,
                                                                     outtree_fname=inferred_full_popn_section_treefile,
                                                                     fastree_exe=FASTTREE_EXE)


        # Calculate HyPhy dN/dS for the full sample_genomes population fasta
        hyphy_filename_prefix = (OUT_SUBS_DIR + os.sep +
                                       FILENAME_PREFIX + ".break.{}_{}".format(section_nuc_start_base1, section_nuc_end_base1))
        hyphy_handler.calc_dnds(codon_fasta_filename=full_popn_section_fasta,
                                tree_filename=inferred_full_popn_section_treefile,
                                hyphy_filename_prefix=hyphy_filename_prefix,
                                threads=PROCS)


    # Combine dN/dS TSV from each recominant section into one TSV
    full_popn_dnds_tsv = OUT_SUBS_DIR + os.sep + FILENAME_PREFIX + ".dnds.tsv"
    with open(full_popn_dnds_tsv, "w") as fh_dnds_out:
        import csv
        writer = csv.DictWriter(fh_dnds_out, fieldnames=hyphy_handler.FIELD_NAMES, delimiter="\t")
        writer.writeheader()

        for i, (section_nuc_start_base1, section_nuc_end_base1) in enumerate(recombo_sections):
            section_dnds_tsv = (OUT_SUBS_DIR + os.sep +
                                FILENAME_PREFIX + ".break.{}_{}.dnds.tsv".format(section_nuc_start_base1, section_nuc_end_base1))
            with open(section_dnds_tsv, 'rU') as fh_in_dnds:
                reader = csv.DictReader(fh_in_dnds, delimiter="\t")
                for inrow in reader:
                    codon_site_wrt_section_0based = int(inrow[hyphy_handler.HYPHY_TSV_SITE])
                    codon_site_wrt_genome_0based = codon_site_wrt_section_0based + ((section_nuc_start_base1 - 1) / 3)

                    outrow = inrow
                    outrow[hyphy_handler.HYPHY_TSV_SITE] = codon_site_wrt_genome_0based
                    writer.writerow(outrow)








