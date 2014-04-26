import os
import errno


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