#pull a slice out of a *.csf file (that our miseq pipeline generates) as a FASTA file.  It works!  And you can build some respectable trees from slices of 300 bp or more.
import sys

def convert_csf (csf_handle):
    """
    Extract the header, offset, and seq from the csf, return a list of
    [header, seq] tuples of decompressed sequence data.
    """
    fasta = []
    for line in csf_handle:
        header, offset, seq = line.strip('\n').split(',')
        fasta.append([header, seq, int(offset), int(offset) + len(seq)])        
    return fasta

try:
    infile = open(sys.argv[1], 'rU')
    fasta = convert_csf(infile)
    infile.close()

    # coordinates defined relative to reference of alignment
    from_left = int(sys.argv[2])
    to_right = int(sys.argv[3])
    span = to_right - from_left

    outfile = open(sys.argv[4], 'w')
except IndexError:
    print 'python extract_csf.py [csf file] [left] [right] [outfile]'
    sys.exit()
except:
    raise

for i, (h, s, left, right) in enumerate(fasta):
    if left > to_right:
        # we are past the window
        # reads are in order, so exit from loop
        break
        
    if right < from_left:
        continue
    
    # pad the read so we are working on a consistent coordinate system
    # and add gaps out to the desired right bound
    s2 = '-'*left + s
    
    if right < to_right:
        # read falls short of desired right bound
        s2 += '-' * (to_right - right)
    
    s3 = s2[from_left:to_right]
    
    overlap = sum(map(lambda x: int(x != '-'), list(s3)))
    
    if overlap / float(span) < 0.75:
        continue
    
    outfile.write('>%d\n%s\n' % (i, s3))

outfile.close()
