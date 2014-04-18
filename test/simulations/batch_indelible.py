"""
Batch INDELible
Write control.txt file and execute indelible to simulate sequences.
"""
import sys
import os

treefile = sys.argv[1]
handle = open(treefile, 'rU')

tree_string = handle.readline()
handle.close()

# factor to stretch entire tree by (mutation rate)
global_scaling_factors = [20.0, 50.0]#[1.0, 2.0, 5.0, 10.0]
#omegas = [0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

# Gamma(shape=1.5, rate=3), histogram with 50 breaks
omegas = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65,
        1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25, 3.35,
        3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15, 4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05,
        5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 6.05]

prop = [0.1038610, 0.1432928, 0.1381003, 0.1212857, 0.1020363, 0.0835798, 0.0673901, 0.0535906,
        0.0422005, 0.0329969, 0.0258732, 0.0200207, 0.0154661, 0.0118681, 0.0090903, 0.0070075,
        0.0053782, 0.0040914, 0.0031212, 0.0023785, 0.0017896, 0.0013684, 0.0010189, 0.0007866,
        0.0005856, 0.0004496, 0.0003366, 0.0002510, 0.0001857, 0.0001455, 0.0001097, 0.0000839,
        0.0000653, 0.0000442, 0.0000391, 0.0000294, 0.0000200, 0.0000160, 0.0000124, 0.0000084,
        0.0000066, 0.0000052, 0.0000031, 0.0000022, 0.0000017, 0.0000016, 0.0000009, 0.0000008,
        0.0000008, 0.0000003, 0.0000002, 0.0000002, 0.0000002, 0.0000001, 0.0000003, 0.0000001,
        0.0000002, 0.0000001, 0.0000001, 0.0000001, 0.0000001]

#prop = 1./len(omegas)
#prop = [0.2468698, 0.03992828, 0.0004993992]
n_sites = 10000 # codons
kappa = 8.0

for scaling_factor in global_scaling_factors:
    handle = open('control.txt', 'w')

    # write minimal contents of INDELible control file
    handle.write('[TYPE] CODON 1\n')
    handle.write('[SETTINGS]\n[output] FASTA\n[printrates] TRUE\n')
    handle.write('[MODEL] M3\n[submodel] %f\n' % kappa)
    
    prop_string = ''
    omega_string = ''
    for i, omega in enumerate(omegas):
        if i < (len(omegas)-1):
            prop_string += ' %f' % prop[i]
        omega_string += ' %1.2f' % omega
    
    handle.write(prop_string + '\n')
    handle.write(omega_string + '\n')
    
    handle.write('[TREE] bigtree %s;\n' % tree_string.rstrip('0123456789.;:\n'))
    handle.write('[treelength] %1.1f\n' % scaling_factor)
    handle.write('[PARTITIONS] partitionname\n')
    handle.write('  [bigtree M3 %d]\n' % n_sites)
    handle.write('[EVOLVE] partitionname 1 scaling_%1.1f\n' % (scaling_factor))

    handle.close()
    
    os.system('indelible')  # indelible
    # miseq - not gemsim

