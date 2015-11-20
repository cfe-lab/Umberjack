#! /usr/bin/env python

import string
import csv
import xml.etree.ElementTree as ET
import sys

from ASGsim import *
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", help="filepath and filename prefix (excluding the file suffix) for the output tree newick file and rate tsv")
parser.add_argument("-n", help="number of leafs", type=int)
parser.add_argument("-r", help="integer seed.  Passed to numpy random generator.  Must be < 2^32", type=int)
parser.add_argument("-s", help="selection rate.  Default=0.01", type=float, default=0.01)
parser.add_argument("-g", help="generations", type=int)
parser.add_argument("-m", help="mutation rate.  Default=0.0001", type=float, default=0.0001)


if len(sys.argv) <= 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

root = args.f
num_leafs = args.n
seed = args.r
mutation_rate = args.m
selection_rate = args.s
generations = args.g

if not root or not num_leafs:
    print "Must specify filepath and filename prefix (-f) and number of leafs (-n)."
    parser.print_help()
    sys.exit(1)

whens = [generations]  # 1 HIV generation = 1 day.  5000 = 13 years
how_manys = [num_leafs]

# first 3 parameters are: resample_rate, selection_rate, mutation_rate
myASG = ASG(0.00001, selection_rate, mutation_rate, whens, how_manys, growth_rate = 0.005)

print("Building ASG...")
myASG.buildASG(random_seed=seed)
print("Assigning types...")
myASG.assign_types()
print("Marking genealogy...")
myASG.mark_true_genealogy()
    
# myASG.Graphviz("GV_test.png")

print("Creating Newick tree... " + os.path.abspath(root+".nwk"))
annot, Newick = myASG.FigTree()

# Write the annotations to FT_test.tsv
f = open(root+".tsv", "wb")
annot_writer = csv.writer(f, delimiter='\t')
annot_writer.writerow(("leaf", "taxa", "type", "mut", "time"))
annot_writer.writerows(annot)
f.close()

# Write the tree out
f = open(root+".nwk", "wb")
f.write(Newick)
f.close()

# Write out the PhyloXML representation
#pxml = myASG.PhyloXML()

#pxml.write("PXML_test.xml")

# Write out the NeXML representation
#nxml = myASG.NeXML()
#nxml.write("NeXML_test.xml")
