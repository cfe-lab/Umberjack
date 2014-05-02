import os
import errno

NUC_PER_CODON = 3

def create_dir_check(new_dir):
    """
    Creates the directory if it doesn't exist.  Raises errors encountered except 'file already exists' error.
    :param str new_dir:  full file path of new directory
    :raises OSError
    """
    try:
        os.makedirs(new_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def get_seq2len(fasta_filename):
    # keep track of reference contig/chromosome and its length
    """
    Map sequence header to sequence length in a dictionary.

    :rtype dict: the seq2len dictionary
    :param str fasta_filename:
    """
    seq2len = {}
    with open(fasta_filename, 'r') as ref_fasta_fh:
        header = ''
        seq_len = 0
        for line in ref_fasta_fh:
            line = line.rstrip().split()[0]  # split by whitespace, only take the first chunk

            if line[0] == '>':
                if seq_len:
                    seq2len[header] = seq_len

                header = line[1:]
                seq_len = 0
            else:
                seq_len += len(line)
        # Reached EOF.  Write out the cached header and sequence length to dictionary.
        seq2len[header] = seq_len
    return seq2len


def get_fasta_headers(fasta_filename):
    """
    Gets a list of headers from the fasta.  Does not include the ">" header prefix.
    Does not include anything after the first whitespace in the header.
    :rtype list[str] : list of headers
    :param fasta_filename : full file path to the fasta file
    """
    headers = []
    with open(fasta_filename, 'r') as fasta_fh:
        for line in fasta_fh:
            if line[0] == '>':
                header = line[1:].rstrip().split()
                headers.extend(header)
    return headers


def get_longest_seq_size_from_fasta(fasta_filename):
    """
    Gets the size of the longest sequence in the fasta.
    :rtype int :  size in bp of the longest sequence in the fasta.  Or -1 if error.
    :param fasta_filename: full filepath to the fasta.
    """
    longest_seq_len = -1
    with open(fasta_filename, 'r') as fasta_fh:
        seq_len = 0
        for line in fasta_fh:
            line = line.rstrip()
            if line[0] == '>':
                longest_seq_len = max(seq_len, longest_seq_len)
                seq_len = 0
            else:
                seq_len += len(line)
        longest_seq_len = max(seq_len, longest_seq_len)

    return longest_seq_len

def get_total_seq_from_fasta(fasta_filename):
    count = 0
    with open(fasta_filename) as fh:
        for line in enumerate(fh):
            if line[0] == '>':
                count += 1
    return count



def get_total_codons_by_pos(msa_fasta_filename):
    """
     Gets List of  total number of sequences with an unambiguous codon for every codon position.
     ASSUME that fasta sequences start at ORF start.
     ASSUME that fasta sequences are multiple sequence aligned

    :rtype : list of codon counts for each position
    :param str msa_fasta_filename:  full filepath to nucleotide fasta
    """
    longest_seq = get_longest_seq_size_from_fasta(msa_fasta_filename)
    total_unambig_codon_by_pos = [0] * ((longest_seq/NUC_PER_CODON) + 1)
    with open(msa_fasta_filename, 'r') as fh:
        seq = ""
        for line in fh:
            line = line.rstrip()
            if line[0] == '>':
                if len(seq) >= 2:
                    for nuc_pos in range(0, len(seq), 3):
                        codon_1st_nuc = seq[nuc_pos].upper()
                        if len(seq) > nuc_pos+1:
                            codon_2nd_nuc = seq[nuc_pos + 1].upper()
                        else:
                            codon_2nd_nuc = "-"
                        if not (codon_1st_nuc == 'N' or codon_1st_nuc == '-' or codon_2nd_nuc == 'N' or codon_2nd_nuc == '-'):
                            codon_pos = nuc_pos/NUC_PER_CODON
                            total_unambig_codon_by_pos[codon_pos] += 1

                seq = ""
            else:
                seq += line

        if len(seq) >= 2:
            for nuc_pos in range(0, len(seq), 3):
                codon_1st_nuc = seq[nuc_pos].upper()
                if len(seq) > nuc_pos+1:
                    codon_2nd_nuc = seq[nuc_pos + 1].upper()
                else:
                    codon_2nd_nuc = "-"
                if not (codon_1st_nuc == 'N' or codon_1st_nuc == '-' or codon_2nd_nuc == 'N' or codon_2nd_nuc == '-'):
                    codon_pos = nuc_pos/NUC_PER_CODON
                    total_unambig_codon_by_pos[codon_pos] += 1

    return total_unambig_codon_by_pos


class Consensus:

    __base_count = {'A':0, 'C':0, 'T':0, 'G':0}
    def __init__(self):
        self.seq = []


    def add_base(self, base, pos_0based):
        base = base.upper()
        if pos_0based >= len(self.seq):
            for i in range(pos_0based - len(self.seq) + 1):
                self.seq.append(Consensus.__base_count.copy())

        if self.seq[pos_0based].has_key(base):
            self.seq[pos_0based][base] += 1

    def get_consensus(self):
        consensus = ""
        for base_count in self.seq:
            max_base = max(base_count, key=base_count.get)
            consensus += max_base
        return consensus

    def print_stats(self):
        for base_count in self.seq:
            print base_count
            if base_count['A'] == base_count['C'] and base_count['G'] == base_count['T'] and base_count['A'] == base_count['T']:
                print " TIE!"
            else:
                print "NARP"
            print "\n"



def get_consensus_from_msa(msa_fasta_filename, consensus_fasta_filename):
    """
    Gets the consensus from a multiple sequence aligned fasta
    """
    consensus = ""
    with open(msa_fasta_filename, 'r') as in_fh, open(consensus_fasta_filename, 'w') as out_fh:
        seq = ""
        consensus = Consensus()
        for line in in_fh:
            line = line.rstrip()
            if line[0] == '>':
                for pos_0based, base in enumerate(seq):
                    consensus.add_base(base, pos_0based=pos_0based)
                seq = ""
            else:
                seq += line

        for pos_0based, base in enumerate(seq):
            consensus.add_base(base, pos_0based=pos_0based)

        consensus_seq = consensus.get_consensus()
        out_fh.write(">consensus " + msa_fasta_filename + "\n")
        out_fh.write(consensus_seq + "\n")

        consensus.print_stats()

