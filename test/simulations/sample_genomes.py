"""
Create genome alignment by random sampling columns from INDELible simulation outputs.
"""
import sys
from seqUtils import *
import random

root = '/Users/apoon/wip/sliders/simulation/'

global_scaling_factors = [10.0, 100.0]

# from batch_indelible.py, to convert from rate class to omega
omegas = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65,
        1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25, 3.35,
        3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15, 4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05,
        5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 6.05]

rate_class_to_omega = dict([(i, omega) for i, omega in enumerate(omegas)])


n_sites = 3000 # target size of genome in codons
#n_taxa = 10000

#selection = random.sample(range(len(fasta)), n_taxa)

# make up a genome with five intervals 
# mutation rate factors 1x, 10x, 100x, 10x, 1x 
# 500, 500, 1000, 500, 500 codons
scaling_factors = [1.0, 2.0, 5.0, 10.0, 100.0, 10.0, 5.0, 2.0, 1.0]
lengths = [100, 200, 200, 500, 1000, 500, 200, 200, 100]  # codons

left = 0  # gets updated with each interval

new_fasta = {}

ratefile = open('sample_genomes.rates', 'w')  # keep track of each codon site omega

for i, scaling_factor in enumerate(scaling_factors):
    # read FASTA
    infile = open(root+'version1/scaling_'+str(scaling_factor)+'.fas', 'rU')
    fasta = convert_fasta(infile.readlines())  # returns a List of header, sequence lists
    infile.close()
    
    # if this is first time, transfer headers
    if len(new_fasta) == 0:
        new_fasta = dict([(h, '') for h, s in fasta])
    
    # transfer intervals
    for h, s in fasta:
        new_fasta[h] += s[ (3*left) : (3*(left+lengths[i])) ]  # in bps
    
    # read rates, starts on line 11 (10 in zero-index)
    infile = open(root+'version1/scaling_'+str(scaling_factor)+'_RATES.txt', 'rU')
    lines = infile.readlines()
    infile.close()
    for line in lines[(10+left):(10+left+lengths[i])]:
        site, rate_class, partition = line.strip('\n').split('\t')
        rate_class = int(rate_class)
        ratefile.write('%d,%1.1f,%d,%f\n' % (i, scaling_factor, rate_class, omegas[rate_class]))
    
    left += lengths[i]  # update left bound
    

# output            
outfile = open('sample_genomes.fas', 'w')

for h in new_fasta.iterkeys():
    outfile.write('>%s\n%s\n' % (h, new_fasta[h]))

outfile.close()


