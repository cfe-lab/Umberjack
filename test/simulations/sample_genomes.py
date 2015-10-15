"""
Create genome alignment by random sampling columns from INDELible simulation outputs.
"""
from collections import OrderedDict
import Utility
import csv
import os
import random
import ctypes
from collections import deque
from argparse import ArgumentParser

BP_PER_CODON = 3

# need to allow random scaling factor per codon
# random assign scaling factor in blocks
# total codon sites in total
parser = ArgumentParser()
parser.add_argument("-scaling_factors", help="Comma delimited list of scaling factors", type=str)
parser.add_argument("-codons_per_block", help="Number of codons per block", type=int)
parser.add_argument("-total_blocks", help="Total blocks in the genome", type=int)
parser.add_argument("-seed", help="random seed for reproducible output", type=int)
parser.add_argument("-input_dir", help="input directory under which INDELible fastas under different mutation rates will exist")
parser.add_argument("-input_filename_prefix", help="filename prefix (no extension or path directory) of the INDELible fastas")
parser.add_argument("-output_dir", help="output directory under which randomized genome will exist")
parser.add_argument("-output_filename_prefix", help="output filename prefix (no extension or path directory) fastas with randomized mutation")

args = parser.parse_args()
scaling_factors = args.scaling_factors.split(",")
output_dir = args.output_dir
output_filename_prefix = args.output_filename_prefix
seed = args.seed
codons_per_block = args.codons_per_block
total_blocks = args.total_blocks
indelible_output_dir = args.input_dir
indelible_filename_prefix = args.input_filename_prefix

if total_blocks % len(scaling_factors) != 0:
    raise ValueError("The number of scaling factors should divide evenly into the number of blocks")


if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ratefile = open(output_dir + os.sep + output_filename_prefix + ".rates.csv", 'w')  # keep track of each codon site omega
ratefile.write('Site,Scaling_factor,Rate_class,Omega\n')

# we randomize which blocks are assigned which mutation rate by shuffling the blocks and
# then assigning mutation rates in order to the shuffled blocks
randomizer = random.Random(seed)
blocks = range(0,total_blocks)
randomizer.shuffle(blocks)
shuffled_blocks_queue = deque(blocks)

new_fasta = None
new_anc_fasta = None
# Iterate through mutation scaling rate in order, but assign it to random block
# Do it this way instead of iterating through blocks and assigning them to random mutation rates,
# so that we limit the number of times we have to hash entire INDELible fastas for each mutation rate.
for scaling_factor in scaling_factors:

    scaling_factor = scaling_factor.lstrip().rstrip()

    # read INDELible tree sequence FASTA
    # <output_filename_prefix>_<scaling_factor>_TRUE.fasta are fastas containing the INDELible tip sequences of a phylogenetic tree
    # The tree mutation rate is scaled by <scaling_factor>.
    infile = open(indelible_output_dir+ os.sep + scaling_factor + os.sep + indelible_filename_prefix + "." +scaling_factor+'_TRUE.fasta', 'rU')
    indelible_fasta = Utility.convert_fasta(infile.readlines())  # returns a List of (header, sequence) tuples
    infile.close()


    # read INDELible ancestor  FASTA
    # <output_filename_prefix>_<scaling_factor>_ANCESTRAL.fasta are fastas containing the INDELible inner node sequences of a phylogenetic tree
    fh_in_anc_fasta_ = open(indelible_output_dir+ os.sep + scaling_factor + os.sep + indelible_filename_prefix + "." +scaling_factor+'_ANCESTRAL.fasta', 'rU')
    indelible_anc_fasta = Utility.convert_fasta(fh_in_anc_fasta_.readlines())  # returns a List of (header, sequence) tuples
    fh_in_anc_fasta_.close()


    # if this is first time, transfer header
    if not new_fasta:
        new_fasta = OrderedDict()
        new_anc_fasta = OrderedDict()
        for h, s in indelible_fasta:
            new_fasta[h] = ctypes.create_string_buffer(s)
        for h, s in indelible_anc_fasta:
            new_anc_fasta[h] = ctypes.create_string_buffer(s)


    # Take codons from the sequences from the INDELible population matching the current mutation scaling rate
    blocks_per_scaling = total_blocks / len(scaling_factors)
    selected_blocks = []
    for i in range(0, blocks_per_scaling):
        # Each item in queue is a 0-based block index wrt original INDELible sequence.
        block_index  = shuffled_blocks_queue.popleft()  # 0-based block index wrt original INDELible sequence
        block_start_bp = BP_PER_CODON * codons_per_block * block_index  # 0-based nucleotide index wrt original INDELible sequence
        block_end_bp = BP_PER_CODON * codons_per_block * (block_index+1)  # 0-based nucleotide index wrt original INDELible sequence
        for h, s in indelible_fasta:
            new_fasta[h][block_start_bp : block_end_bp] = s[block_start_bp : block_end_bp]
        for h, s in indelible_anc_fasta:
            new_anc_fasta[h][block_start_bp : block_end_bp] = s[block_start_bp : block_end_bp]

        selected_blocks.append(block_index)

    indelible_rates_csv = indelible_output_dir + os.sep + scaling_factor + os.sep +  indelible_filename_prefix + "." + scaling_factor+'_RATES.csv'
    with open(indelible_rates_csv, 'rU') as fh_rates_in:
        reader = csv.DictReader(fh_rates_in)  # Columns:  Site	Class	Partition	Inserted?	Omega
        for line in reader:
            codonsite_1based = int(line["Site"])
            orig_block = (codonsite_1based - 1) / codons_per_block
            if orig_block in selected_blocks:
                ratefile.write('%d,%s,%s,%s\n' % (codonsite_1based, scaling_factor, line["Class"], line["Omega"]))



# output
out_fasta_fname = output_dir+ os.sep + output_filename_prefix + ".fasta"
out_anc_fasta_fname = output_dir+ os.sep + output_filename_prefix + ".ancestral.fasta"
with open(out_fasta_fname, 'w') as fh_out_fasta:
    for h in new_fasta.iterkeys():
        fh_out_fasta.write('>%s\n%s\n' % (h, new_fasta[h].value.rstrip()))  # Need rstrip since extra whitespace added to .value

with open(out_anc_fasta_fname, 'w') as fh_out_anc_fasta:
    for h in new_anc_fasta.iterkeys():
        fh_out_anc_fasta.write('>%s\n%s\n' % (h, new_anc_fasta[h].value.rstrip()))  # Need rstrip since extra whitespace added to .value



# output consensus
Utility.write_consensus_from_msa(out_fasta_fname, out_fasta_fname.replace(".fasta", ".consensus.fasta"))


