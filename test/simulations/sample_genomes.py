"""
Create genome alignment by random sampling columns from INDELible simulation outputs.
"""
import collections
import Utility
import csv
import sys
import os

BASES_PER_CODON = 3

root = sys.argv[1]
input_prefix = sys.argv[2]  # filename prefix (no parent directory) of the indelible files before the mutation rate
output_prefix = sys.argv[3]


Interval = collections.namedtuple(typename='Interval', field_names=['scaling_factor', 'num_codons'], verbose=True)
interval_csv = root + os.sep + "intervals.csv"
intervals = []
with open(interval_csv, 'rU') as fh_intervals:
    reader = csv.DictReader(fh_intervals)
    for row in reader:
        intervals.extend([Interval(scaling_factor=float(row["scaling_factor"]), num_codons=int(row["num_codons"]))])

interval_start_codon_1based = 1 # gets updated with each interval
new_codon_site_1based = 1
new_fasta = {}


ratefile = open(root + os.sep + output_prefix + ".rates.csv", 'w')  # keep track of each codon site omega
ratefile.write('Site,Interval,Scaling_factor,Rate_class,Omega\n')

for interval_idx, interval in enumerate(intervals):
    scaling_factor = interval.scaling_factor
    num_codons = interval.num_codons
    # read INDELible tree sequence FASTA
    # scaling_<scaling_factor>_TRUE.fas are fastas containing the INDELible node sequences of a phylogenetic tree
    # The tree mutation rate is scaled by <scaling_factor>.
    infile = open(root+ os.sep + input_prefix + str(scaling_factor)+'_TRUE.fasta', 'rU')
    fasta = Utility.convert_fasta(infile.readlines())  # returns a List of (header, sequence) tuples
    infile.close()
    
    # if this is first time, transfer header
    if len(new_fasta) == 0:
        new_fasta = dict([(h, '') for h, s in fasta])
    
    # transfer intervals
    for h, s in fasta:
        new_fasta[h] += s[ BASES_PER_CODON*(interval_start_codon_1based-1) : BASES_PER_CODON*(interval_start_codon_1based+num_codons-1) ]  # in bps

    indelible_rates_csv = root + os.sep +  input_prefix + str(scaling_factor)+'_RATES.csv'
    with open(indelible_rates_csv, 'rU') as fh_rates_in:
        reader = csv.DictReader(fh_rates_in)  # Columns:  Site	Class	Partition	Inserted?	Omega
        for line in reader:
            if int(line["Site"]) < interval_start_codon_1based:
                continue
            if int(line["Site"]) >= interval_start_codon_1based + num_codons:
                interval_start_codon_1based = interval_start_codon_1based + num_codons
                break
            ratefile.write('%d,%d,%1.1f,%s,%s\n' % (new_codon_site_1based, interval_idx, scaling_factor, line["Class"], line["Omega"]))
            new_codon_site_1based += 1


# output
out_fasta_fname = root+ os.sep + output_prefix + ".fasta"
out_fasta_file = open(out_fasta_fname, 'w')



for h in new_fasta.iterkeys():
    out_fasta_file.write('>%s\n%s\n' % (h, new_fasta[h]))

out_fasta_file.close()

# output consensus
Utility.get_consensus_from_msa(out_fasta_fname, out_fasta_fname.replace(".fasta", ".consensus.fasta"))


