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

    sections = []
    if num_breaks == 0:
        sectn_start_nuc_base1 = 1
        sectn_end_nuc_base1 = genome_codons * 3
        sections.append((sectn_start_nuc_base1, sectn_end_nuc_base1))

    # We define break_site as the codon position immediately before the strand switch in recombination
    # last position can't be breakpoint because there is no position immediately after it to switch strands to
    # 0-based positions
    breakpt_codons = sorted(randomizer.sample(range(1, genome_codons), num_breaks))

    for b, breakpt_codon_base0 in enumerate(breakpt_codons):
        # Although breakpoints are given in 0-based codon positions,
        # we follow the Umbjerack custom of specifying genome partitions within filenames with 1-based nucleotide positions
        breakpt_nuc_base1 = (breakpt_codon_base0 * 3) + 1
        if b == 0:  # first breakpoint
            sectn_start_nuc_base1 = 1
        else:
            sectn_start_nuc_base1 = breakpt_nuc_base1

        if b == len(breakpt_codons) - 1:  # last breakpoint
            sectn_end_nuc_base1 = genome_codons * 3
        else:
            next_breakpt_codon_base0 = breakpt_codons[b+1]
            next_breakpoint_nuc_base1 = (next_breakpt_codon_base0 * 3) + 1
            sectn_end_nuc_base1 = next_breakpoint_nuc_base1 - 1

        sections.append((sectn_start_nuc_base1,  sectn_end_nuc_base1))

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
