"""
Prune terminal branches on phylogeny given minimum distance cutoff.
"""
import sys

from Bio import Phylo

from test.simulations.seqUtils import convert_fasta


try:
    fasta = convert_fasta(open(sys.argv[1], 'rU').readlines())
    t = Phylo.read(sys.argv[2], 'newick')
    cutoff = float(sys.argv[3])
    outfile = open(sys.argv[4], 'w')
except:
    print 'Prune terminal branches on phylogeny given minimum distance cutoff.'
    print 'Filter the original sequence alignment based on the pruned tree.'
    print 'python prune_newick.py [fasta] [newick] [cutoff] [outfile]'
    raise


# make sure that tip names match sequence names in FASTA
seq_names = [h for h, s in fasta]
seq_names.sort()
tip_names = [tip.name for tip in t.get_terminals()]
tip_names.sort()

assert tip_names == seq_names, 'tip names do not match'
# actually this won't work because tip names got truncated :-P

while True:
    tips = t.get_terminals()
    pruned = False
    for tip in tips:
        if tip.branch_length < cutoff:
            #print 'pruning', tip.name
            t.prune(tip)
            # reload tree
            pruned = True
            break
    if not pruned:
        break

#Phylo.write(t, file=sys.argv[1]+'.pruned', format='newick')
print 'pruned tree to', len(tips), 'tips'
keep_tips = [tip.name for tip in tips]
for h, s in fasta:
    if h in keep_tips:
        outfile.write('>%s\n%s\n' % (h, s))

outfile.close()

