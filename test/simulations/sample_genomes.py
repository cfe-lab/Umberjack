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

interval_csv = sys.argv[1]
output_dir = sys.argv[2]
output_prefix = sys.argv[3]
seed = int(sys.argv[4])



new_fasta = {}

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ratefile = open(output_dir + os.sep + output_prefix + ".rates.csv", 'w')  # keep track of each codon site omega
ratefile.write('Site,Interval,Scaling_factor,Rate_class,Omega\n')


# Count number of codons all together in all intervals
total_codon_sites = 0
with open(interval_csv, 'rU') as fh_intervals:  #scaling_factor,num_codons,outdir,out_filename_prefix
    reader = csv.DictReader(fh_intervals)
    for interval_row in reader:
        total_codon_sites += int(interval_row["num_codons"])

randomizer = random.Random(seed)
random_codons = range(0,total_codon_sites)
randomizer.shuffle(random_codons)  # we randomize which codons are assigned which mutation rate
random_codons_queue = deque(random_codons)

#for interval_idx, interval in enumerate(intervals):
with open(interval_csv, 'rU') as fh_intervals:
    reader = csv.DictReader(fh_intervals)
    for interval_idx, row in enumerate(reader):


        scaling_factor = float(row["scaling_factor"])
        interval_num_codons = int(row["num_codons"])
        indelible_output_dir = row["outdir"]
        indelible_output_prefix = row["out_filename_prefix"]



        # read INDELible tree sequence FASTA
        # scaling_<scaling_factor>_TRUE.fas are fastas containing the INDELible node sequences of a phylogenetic tree
        # The tree mutation rate is scaled by <scaling_factor>.
        infile = open(indelible_output_dir+ os.sep + indelible_output_prefix + str(scaling_factor)+'_TRUE.fasta', 'rU')
        indelible_fasta = Utility.convert_fasta(infile.readlines())  # returns a List of (header, sequence) tuples
        infile.close()

        # if this is first time, transfer header
        if len(new_fasta) == 0:
            new_fasta = OrderedDict()
            for h, s in indelible_fasta:
                new_fasta[h] = ctypes.create_string_buffer(s)


        # transfer intervals
        interval_random_codon_sites_1based = []
        for i in range(0, interval_num_codons):
            random_codon_site  = random_codons_queue.popleft()  #0-based random codon sites for this mutation rate
            interval_random_codon_sites_1based.append(random_codon_site+1)
            for h, s in indelible_fasta:
                new_fasta[h][BASES_PER_CODON*random_codon_site : BASES_PER_CODON*(random_codon_site+1)] = s[BASES_PER_CODON*random_codon_site : BASES_PER_CODON*(random_codon_site+1)]

        indelible_rates_csv = indelible_output_dir + os.sep +  indelible_output_prefix + str(scaling_factor)+'_RATES.csv'
        with open(indelible_rates_csv, 'rU') as fh_rates_in:
            reader = csv.DictReader(fh_rates_in)  # Columns:  Site	Class	Partition	Inserted?	Omega
            for line in reader:
                site_1based = int(line["Site"])
                if site_1based in interval_random_codon_sites_1based:
                    ratefile.write('%d,%d,%1.1f,%s,%s\n' % (site_1based, interval_idx, scaling_factor, line["Class"], line["Omega"]))



# output
out_fasta_fname = output_dir+ os.sep + output_prefix + ".fasta"
out_fasta_file = open(out_fasta_fname, 'w')



for h in new_fasta.iterkeys():
    out_fasta_file.write('>%s\n%s\n' % (h, new_fasta[h].value.rstrip()))  # Need rstrip since extra whitespace added to .value

out_fasta_file.close()

# output consensus
Utility.get_consensus_from_msa(out_fasta_fname, out_fasta_fname.replace(".fasta", ".consensus.fasta"))


