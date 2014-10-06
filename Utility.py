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
    """
    Keep track of reference contig/chromosome and its length
    :return: the sequence to length dictionary
    :rtype : dict {str:int}
    :param str fasta_filename: full filepath to fasta
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
    :return: list of headers
    :rtype:  list[str]
    :param str fasta_filename : full file path to the fasta file
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
    :return: size in bp of the longest sequence in the fasta.  Or -1 if error.
    :rtype : int
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
    """
    Gets the total sequences from the fasta file.
    :return:  total sequences in the fasta file.
    :rtype : int
    :param fasta_filename: full filepath to fasta
    """
    count = 0
    with open(fasta_filename, 'rU') as fh:
        for line in fh:
            if line[0] == '>':
                count += 1
    return count

def get_total_nongap_nuc_by_pos(msa_fasta_filename, pos):
    """
    Gets list where each element is the total number of sequences with nongap at that position.
    ASSUME that fasta sequences start at ORF start.
    ASSUME that fasta sequences are multiple sequence aligned.

    :param msa_fasta_filename:
    :param int pos:  0-based position index into multiple sequence alignment columns.
    :return:
    """
    pos_total_nongap = 0
    with open(msa_fasta_filename, 'r') as fh:
        seq = ""
        for line in fh:
            line = line.rstrip()
            if line[0] == '>' and seq:
                pos_total_nongap += 1 if seq[pos] != "N" and seq[pos] != "n" and seq[pos] != "-" else 0
                seq = ""
            else:
                seq += line

        pos_total_nongap += 1 if seq[pos] != "N" and seq[pos] != "n" and seq[pos] != "-" else 0

    return pos_total_nongap

def get_total_nongap_nuc_all_pos(msa_fasta_filename):
    """
    Gets list where each element is the total number of sequences with nongap at that position.
    ASSUME that fasta sequences start at ORF start.
    ASSUME that fasta sequences are multiple sequence aligned.

    :param msa_fasta_filename:
    :return:
    """

    longest_seq_size = get_longest_seq_size_from_fasta(msa_fasta_filename)
    total_nongap_by_pos = [0] * longest_seq_size
    with open(msa_fasta_filename, 'r') as fh:
        seq = ""
        for line in fh:
            line = line.rstrip()
            if line[0] == '>':
                for pos, nuc in enumerate(seq):
                    total_nongap_by_pos[pos] += 1 if nuc != "N" and nuc != "n" and nuc != "-" else 0

                seq = ""
            else:
                seq += line

        for pos, nuc in enumerate(seq):
            total_nongap_by_pos[pos] += 1 if nuc != "N" and nuc != "n" and nuc != "-" else 0

    return total_nongap_by_pos


def get_total_codons_by_pos(msa_fasta_filename):
    """
     Gets List of  total number of sequences with an unambiguous codon for every codon position.
     ASSUME that fasta sequences start at ORF start.
     ASSUME that fasta sequences are multiple sequence aligned

    :return: list of codon counts for each position
    :rtype : list of int
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


class _Consensus:
    """
    Internal use only.  Keeps track of consensus information for a sequence.
    """

    __base_count = {'A': 0, 'C': 0, 'T': 0, 'G': 0}

    def __init__(self):
        self.seq = []


    def add_base(self, base, pos_0based):
        """
        Add a base at the given nucleotide position.

        :param str base:  nucleotide base
        :param int pos_0based: 0-based nucleotide position
        """
        base = base.upper()
        if pos_0based >= len(self.seq):
            for i in range(pos_0based - len(self.seq) + 1):
                self.seq.append(_Consensus.__base_count.copy())

        if self.seq[pos_0based].has_key(base):
            self.seq[pos_0based][base] += 1

    def get_consensus(self):
        """
        Return the consensus for the entire sequence
        :rtype : str
        """
        consensus = ""
        for base_count in self.seq:
            max_base = max(base_count, key=base_count.get)
            consensus += max_base
        return consensus

    def print_stats(self):
        """
        Print information stats to stdout about the set of sequences.
        """
        total_mut = 0
        total_bases = 0
        for nucpos, base_count in enumerate(self.seq, start=1):
            mut = sum(base_count.values()) - max(base_count.values())
            total_bases += sum(base_count.values())
            total_mut += mut
            print "1basedNucPos=" + str(nucpos) + ", mut=" + str(mut) + ", basecounts=",
            print base_count,
            if base_count['A'] == base_count['C'] and base_count['G'] == base_count['T'] and base_count['A'] == base_count['T']:
                print " MAX-TIE!"
            else:
                print " ONEMAX"
        print "Ave Mutations per sequence = " + str(float(total_mut)/len(self.seq))
        print "Ave Mutations per base per sequence = " + str(float(total_mut)/total_bases)


def get_consensus_from_msa(msa_fasta_filename, consensus_fasta_filename):
    """
    Gets the consensus from a multiple sequence aligned fasta and prints out stats to stdout

    :param msa_fasta_filename: full filepath to multiple sequence aligned fasta
    :param consensus_fasta_filename:  full filepath to censusus fasta output
    """
    consensus = ""
    with open(msa_fasta_filename, 'r') as in_fh, open(consensus_fasta_filename, 'w') as out_fh:
        seq = ""
        consensus = _Consensus()
        for line in in_fh:
            line = line.rstrip()
            if line:
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


def convert_fasta (lines):
    """
    Return contents of fasta file as list of tuples (header, sequence)
    :rtype : list of tuples [[str, str]]
    :param lines: list of lines in fasta file
    """
    blocks = []
    sequence = ''
    for i in lines:
        if i[0] == '$': # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                blocks.append([h,sequence])
                sequence = ''	# reset containers
                h = i.strip('\n')[1:]
            else:
                h = i.strip('\n')[1:]
        else:
            sequence += i.strip('\n')
    try:
        blocks.append([h,sequence])	# handle last entry
    except:
        print lines
        raise
    return blocks

# # Only handles the start positions for now
# def convert_1nuc_to_1codon(nuc_startpos_1based):
#     return nuc_startpos_1based/NUC_PER_CODON + 1
#
# # Only handles the start positions for now
# def convert_1codon_to_1nuc(codon_startpos_1based):
#     return ((codon_startpos_1based - 1) * NUC_PER_CODON) + 1