"""
Create genome alignment by random sampling columns from INDELible simulation outputs.
"""
from collections import OrderedDict
import Utility
import csv
import sys
import os
import random
import ctypes
from collections import deque

BASES_PER_CODON = 3

scaling_factors = sys.argv[1].split(",")
output_dir = sys.argv[2]
output_filename_prefix = sys.argv[3]
seed = int(sys.argv[4])
total_codon_sites = int(sys.argv[5])
indelible_output_dir = sys.argv[6]
indelible_filename_prefix = sys.argv[7]





if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ratefile = open(output_dir + os.sep + output_filename_prefix + ".rates.csv", 'w')  # keep track of each codon site omega
ratefile.write('Site,Scaling_factor,Rate_class,Omega\n')


randomizer = random.Random(seed)
random_codons = range(0,total_codon_sites)
randomizer.shuffle(random_codons)  # we randomize which codons are assigned which mutation rate
random_codons_queue = deque(random_codons)

new_fasta = {}
for scaling_factor in scaling_factors:
    num_codons_per_scaling = total_codon_sites / len(scaling_factors)

    scaling_factor = scaling_factor.lstrip().rstrip()

    # read INDELible tree sequence FASTA
    # <output_filename_prefix>_<scaling_factor>_TRUE.fasta are fastas containing the INDELible tip sequences of a phylogenetic tree
    # The tree mutation rate is scaled by <scaling_factor>.
    infile = open(indelible_output_dir+ os.sep + scaling_factor + os.sep + indelible_filename_prefix + "." +scaling_factor+'_TRUE.fasta', 'rU')
    indelible_fasta = Utility.convert_fasta(infile.readlines())  # returns a List of (header, sequence) tuples
    infile.close()

    # if this is first time, transfer header
    if len(new_fasta) == 0:
        new_fasta = OrderedDict()
        for h, s in indelible_fasta:
            new_fasta[h] = ctypes.create_string_buffer(s)


    # Randomly select mutation scaling rate for each codon site
    # Take codons from the sequences from the INDELible population matching that mutation scaling rate
    random_codon_sites_1based = []
    for i in range(0, num_codons_per_scaling):
        random_codon_site  = random_codons_queue.popleft()  #0-based random codon sites for this mutation rate
        random_codon_sites_1based.append(random_codon_site+1)
        for h, s in indelible_fasta:
            new_fasta[h][BASES_PER_CODON*random_codon_site : BASES_PER_CODON*(random_codon_site+1)] = s[BASES_PER_CODON*random_codon_site : BASES_PER_CODON*(random_codon_site+1)]

    indelible_rates_csv = indelible_output_dir + os.sep + scaling_factor + os.sep +  indelible_filename_prefix + "." + scaling_factor+'_RATES.csv'
    with open(indelible_rates_csv, 'rU') as fh_rates_in:
        reader = csv.DictReader(fh_rates_in)  # Columns:  Site	Class	Partition	Inserted?	Omega
        for line in reader:
            site_1based = int(line["Site"])
            if site_1based in random_codon_sites_1based:
                ratefile.write('%d,%s,%s,%s\n' % (site_1based, scaling_factor, line["Class"], line["Omega"]))



# output
out_fasta_fname = output_dir+ os.sep + output_filename_prefix + ".fasta"
with open(out_fasta_fname, 'w') as fh_out_fasta:
    for h in new_fasta.iterkeys():
        fh_out_fasta.write('>%s\n%s\n' % (h, new_fasta[h].value.rstrip()))  # Need rstrip since extra whitespace added to .value



# output consensus
Utility.write_consensus_from_msa(out_fasta_fname, out_fasta_fname.replace(".fasta", ".consensus.fasta"))


