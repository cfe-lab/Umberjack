import os
import errno
import math
import re
import Bio.Phylo  as Phylo

NUC_PER_CODON = 3
STOP_AA = "*"
CODON2AA = {
    'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
    'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
    'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
    'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

    'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
    'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
    'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
    'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

    'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
    'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
    'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
    'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

    'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
    'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
    'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
    'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G',

    # codons with mixture codes
    'GCN' : 'A',
    'TTY' : 'F',
    'AAY' : 'N',
    'GAY' : 'D',
    'TGY' : 'C',
    'CAR' : 'Q',
    'GAR' : 'E',
    'GGN' : 'G',
    'CAY' : 'H',
    'ATH' : 'I',
    'AAR' : 'K',
    'YTR' : 'L','CTN' : 'L',
    'CCN' : 'P',
    'CGN' : 'R','MGR' : 'R',
    'TCN' : 'S','AGY' : 'S',
    'ACN' : 'T',
    'TAY' : 'Y',
    'GTN' : 'V',
    'TAR' : '*','TRA' : '*',
    }

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


def get_total_nongap_all_pos(msa_fasta_filename, gap_chars=["N", "n", "-"]):
    """
    Gets the total sequences that do not contain a gap character at each position
    :param msa_fasta_filename: multiple sequence aligned fasta
    :param gap_chars: list of gap characters
    :return list:  list of total sequences without gap character at each position.  Each element in the list is a position in the msa fasta.
    """
    longest_seq_size = get_longest_seq_size_from_fasta(msa_fasta_filename)
    total_nongap_by_pos = [0] * longest_seq_size
    with open(msa_fasta_filename, 'r') as fh:
        seq = ""
        for line in fh:
            line = line.rstrip()
            if line[0] == '>':
                for pos, character in enumerate(seq):
                    total_nongap_by_pos[pos] += 1 if character not in  gap_chars else 0

                seq = ""
            else:
                seq += line

        for pos, character in enumerate(seq):
            total_nongap_by_pos[pos] += 1 if character not in  gap_chars else 0

    return total_nongap_by_pos


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
    total_unambig_codon_by_pos = [0] * int(math.floor(float(longest_seq)/NUC_PER_CODON)) # if the last codon doesn't have enuf chars, then hyphy ignores it
    with open(msa_fasta_filename, 'rU') as fh:
        seq = ""
        for line in fh:
            line = line.rstrip()
            if line[0] == '>':
                seq = seq.upper().replace("-", "N")
                for nuc_pos in range(0, len(seq), 3):
                    codon = seq[nuc_pos:nuc_pos + NUC_PER_CODON]
                    #codon += ("N" * (NUC_PER_CODON - len(codon)))  # right pad with N's
                    if len(codon) == NUC_PER_CODON and CODON2AA.get(codon):
                        codon_pos = nuc_pos/NUC_PER_CODON
                        total_unambig_codon_by_pos[codon_pos] += 1
                seq = ""
            else:
                seq += line

        seq = seq.upper().replace("-", "N")
        for nuc_pos in range(0, len(seq), 3):
            codon = seq[nuc_pos:nuc_pos + NUC_PER_CODON]
            #codon += ("N" * (NUC_PER_CODON - len(codon)))  # right pad with N's
            if len(codon) == NUC_PER_CODON and CODON2AA.get(codon):
                codon_pos = nuc_pos/NUC_PER_CODON
                total_unambig_codon_by_pos[codon_pos] += 1

    return total_unambig_codon_by_pos


class Consensus:
    """
    Keeps track of consensus information for a sequence.
    """

    __base_count = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N':0, '-':0, 'X':0}

    # ? means unkknown amino acid.  X means left or right pad gap in amino acid
    __aa_count = {'?':0, 'A':0, 'B':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0,
                     'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0,
                     'T':0, 'V':0, 'W':0, 'Y':0, '*':0, 'X':0}
    TRUE_BASES = ["A", "C", "G", "T"]
    AMBIG_BASES = ["N"]
    GAPS = ["-"]
    PADS = ["X"]  # left or right pad
    NON_BASES = AMBIG_BASES + GAPS + PADS


    def __init__(self):
        self.seq = []
        self.aa_seq = []  # Assumes that the sequence starts on ORF

    def parse(self, msa_fasta_filename):
        """
        Reads in multiple sequence alignment file.
        """
        with open(msa_fasta_filename, 'r') as in_fh:
            seq = ""
            for line in in_fh:
                line = line.rstrip()
                if line:
                    if line[0] == '>':
                        if seq:
                            self.add_seq(seq)
                            seq = ""
                    else:
                        seq += line
            if seq:
                self.add_seq(seq)


    def add_seq(self, seq):
        """
        Add a sequence to the consensus
        """

        # Collect the Nucleotide level stats

        # Find the position of the first non-gap char.  Anything before this is a left-pad gap.
        truebase_start = re.search(r"[^\-]", seq).start()
        # Find the position of the last non-gap char.  Anything after this is a right-pad gap.
        truebase_end = re.search(r"[^\-][\-]*$", seq).start()

        for pos_0based, base in enumerate(seq):
            if truebase_start <= pos_0based <= truebase_end:
                self.add_base(base, pos_0based=pos_0based)  # Inner gap, or ACGT, or N
            else:
                self.add_base("X", pos_0based=pos_0based)  # left or right pad gap

        # # Collect the codon level stats
        # truecodon_start = truebase_start/NUC_PER_CODON
        # truecodon_end = truebase_end/NUC_PER_CODON
        #
        # for nuc_pos in range(0, len(seq), 3):
        #     codon = seq[nuc_pos:nuc_pos + NUC_PER_CODON]
        #     codon += ("N" * (NUC_PER_CODON - len(codon)))  # right pad with N's
        #     codon_pos = nuc_pos/NUC_PER_CODON
        #
        #     if truecodon_start <= codon_pos <= truecodon_end:
        #         self.add_codon(codon, codon_pos_0based=codon_pos)   # actual codon
        #     else:
        #         self.add_codon("X", codon_pos_0based=codon_pos)  # left or right pad gap


    def add_codon(self, codon, codon_pos_0based):
        """
        :param codon:
        :param codon_pos_0based:  0-based codon position
        :return:
        """
        # Replace the gaps,  pads with N's so that we can translate the codons.
        codon = codon.upper().replace("-", "N")

        if codon_pos_0based >= len(self.aa_seq):
            for i in range(codon_pos_0based - len(self.aa_seq) + 1):
                self.seq.append(Consensus.__base_count.copy())

        aa = CODON2AA.get(codon, "?")  # TODO:  don't hardcode this.  ? means unknown aa

        if self.aa_seq[codon_pos_0based].has_key(aa):
            self.aa_seq[codon_pos_0based][aa] += 1


    def add_base(self, base, pos_0based):
        """
        Add a base at the given nucleotide position.

        :param str base:  nucleotide base
        :param int pos_0based: 0-based nucleotide position
        """
        base = base.upper()
        if pos_0based >= len(self.seq):
            for i in range(pos_0based - len(self.seq) + 1):
                self.seq.append(Consensus.__base_count.copy())

        if self.seq[pos_0based].has_key(base):
            self.seq[pos_0based][base] += 1


    def get_consensus(self):
        """
        Return the consensus for the entire sequence
        :rtype : str
        """
        consensus = ""
        for base_count in self.seq:
            max_base_count = 0
            max_base = None
            for base, count in base_count.iteritems():
                if max_base_count < count and base not in Consensus.NON_BASES:
                    max_base = base
                    max_base_count = count
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


    def get_conserve(self, pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        :param int pos_0based: 0-based position in the multiple sequence alignment
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return: the fraction of conserved sequences at every position in the multiple sequence alignment.
        :rtype float
        """
        conserve_count = 0
        total_count = 0

        for base, count in self.seq[pos_0based].iteritems():
            if base not in Consensus.NON_BASES:
                total_count += count

                if is_count_ambig:
                    count += self.seq[pos_0based]["N"]/4.0
                if is_count_gaps:
                    count += self.seq[pos_0based]["-"]/4.0
                if is_count_pad:
                    count += self.seq[pos_0based]["X"]/4.0

                if conserve_count < count:
                    conserve_count = count

            elif is_count_ambig and base == "N":
                total_count += count
            elif is_count_gaps and base == "-":
                total_count += count
            elif is_count_pad and base == "X":
                total_count += count

        if total_count:
            return float(conserve_count)/total_count
        else:
            return None


    def get_ave_conserve(self, start_pos_0based, after_end_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        Gets the average per-site conservation across the given range.
        :param start_pos_0based:  0-based nucleotide start position
        :param after_end_pos_0based:  0-based nucleotide end position + 1
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return:  average per-site fraction of conservation across the range.
        """
        total_conserve = 0.0
        for pos in range(start_pos_0based, after_end_pos_0based):
            conserve = self.get_conserve(pos, is_count_ambig, is_count_gaps, is_count_pad)
            total_conserve += conserve if conserve else 0

        ave_conserve = total_conserve / (after_end_pos_0based - start_pos_0based)
        return ave_conserve



    def get_shannon_entropy(self, pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        Gets the Shannon Entropy (measure of bits required to represent each symbol) at the given site.  Only considers A, C, G, T.
        If there are N or -, then adds 1 to each A, C, G, T count.
        :param int pos_0based: 0-based nucleotide position in the multiple sequence alignment
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return: Shannon Entropy.  Log2 based.
        :rtype float
        """
        total_seqs = sum([count for letter, count in self.seq[pos_0based].iteritems()
                          if (is_count_ambig or letter != "N") and (is_count_gaps or letter != "-") and (is_count_pad or letter != "X")])
        if not total_seqs:
            return None

        total_entropy = 0.0
        for letter, count in  self.seq[pos_0based].iteritems():
            if letter not in self.NON_BASES:
                if is_count_ambig:
                    count += self.seq[pos_0based]["N"]/4.0
                if is_count_gaps:
                    count += self.seq[pos_0based]["-"]/4.0
                if is_count_pad:
                    count += self.seq[pos_0based]["X"]/4.0

                if count:
                    p_letter = count / float(total_seqs)  # probability of this letter occuring at this position
                    log_p_letter = math.log(p_letter, 2)  # Log2  probability of letter
                    total_entropy += (p_letter * log_p_letter)

        total_entropy = -total_entropy
        return total_entropy

    def get_ave_shannon_entropy(self, start_pos_0based, after_end_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        Gets the average Shannon Entropy (measure of bits required to represent each symbol) at the given range.
        :param start_pos_0based:  0-based nucleotide start position
        :param after_end_pos_0based:  0-based nucleotide end position + 1
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return: Shannon Entropy.  Log2 based.
        :rtype float
        """
        total_entropy = 0.0
        for pos in range(start_pos_0based, after_end_pos_0based):
            entropy = self.get_shannon_entropy(pos, is_count_ambig, is_count_gaps, is_count_pad)
            total_entropy += entropy if entropy else 0

        ave_entropy = total_entropy / (after_end_pos_0based - start_pos_0based)
        return ave_entropy


    def get_ave_metric_entropy(self, start_pos_0based, after_end_pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        Gets the average per-site metric entropy across the given range.
        :param start_pos_0based:  0-based nucleotide start position
        :param after_end_pos_0based:  0-based nucleotide end position + 1
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return:  average per-site metric entropy across the range.
        """
        total_entropy = 0.0
        for pos in range(start_pos_0based, after_end_pos_0based):
            entropy = self.get_metric_entropy(pos, is_count_ambig, is_count_gaps, is_count_pad)
            total_entropy += entropy if entropy else 0

        ave_entropy = total_entropy / (after_end_pos_0based - start_pos_0based)
        return ave_entropy


    def get_metric_entropy(self, pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        Gets the Metric Shannon Entropy at the given site.  Only considers A, C, G, T.
        :param int pos_0based: 0-based position in the multiple sequence alignment
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return: Metric Shannon Entropy  (the Shannon Entropy divided by log(the number of sequences), \
                which can be compared across sites to measure randomness
        :rtype float
        """
        total_seqs = sum([count for letter, count in self.seq[pos_0based].iteritems()
                          if (is_count_ambig or letter != "N") and (is_count_gaps or letter != "-") and (is_count_pad or letter != "X")])
        if not total_seqs:
            return None
        shannon_entropy = self.get_shannon_entropy(pos_0based, is_count_ambig, is_count_gaps, is_count_pad)
        metric_entropy = shannon_entropy/math.log(total_seqs)
        return metric_entropy


    def get_alignment_len(self):
        """
        :returns: the length of the longest sequence in the alignment
        :rtype: int
        """
        return len(self.seq)

    def get_depth(self, pos_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False):
        """
        :param int pos_0based: 0-based nucleotide position in the multiple sequence alignment
        :param bool is_count_ambig:  whether to include N as 0.25 of A, C, G, T
        :param bool is_count_gaps:  whether to include inner "-" as 0.25 of A, C, G, T
        :param bool is_count_pad:  whether to include outer - as 0.25 of A, C, G, T
        :return: the total sequences with a valid nucleotide (A, C, G, T) at the given position.
        :rtype: int
        """
        total_seqs = sum([count for letter, count in self.seq[pos_0based].iteritems()
                          if (is_count_ambig or letter != "N") and (is_count_gaps or letter != "-") and (is_count_pad or letter != "X")])
        return total_seqs


    def get_ambig_count(self, pos_0based):
        return self.seq[pos_0based]["N"]  # TODO:  don't hardcode

    def get_pad_count(self, pos_0based):
        """
        Here, pad refers to gaps external to the sequence.  I.e.  left and right pad gaps.
        """
        return self.seq[pos_0based]["X"]  # TODO:  don't hardcode

    def get_gap_count(self, pos_0based):
        """
        Here, gap refers to gap internal in the sequence.  I.e.  surrounded by N, A, C, G, T on either side.
        """
        return self.seq[pos_0based]["X"]  # TODO:  don't hardcode

    def get_codon_depth(self, codon_pos_0based, is_count_ambig=False, is_count_pad=False):
        total_seqs = sum([count for letter, count in self.codon_seq[codon_pos_0based].iteritems()
                          if (is_count_ambig or letter != "?") and (is_count_pad or letter != "X")])
        return total_seqs



def get_consensus_from_msa(msa_fasta_filename, consensus_fasta_filename):
    """
    Gets the consensus from a multiple sequence aligned fasta and prints out stats to stdout

    :param msa_fasta_filename: full filepath to multiple sequence aligned fasta
    :param consensus_fasta_filename:  full filepath to censusus fasta output
    """
    consensus = ""
    with open(msa_fasta_filename, 'r') as in_fh, open(consensus_fasta_filename, 'w') as out_fh:
        seq = ""
        consensus = Consensus()
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

        #consensus.print_stats()



def convert_fasta (lines):
    """
    Return contents of fasta file as list of tuples (header, sequence)
    :rtype : list of tuples [(str, str)]
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


def get_seq_dict(fasta):
    """
    :param str fasta:  path to fasta
    :return  dict  {str: str} :  sequence dict  {sequence name: sequence value}
    """
    header2seq = dict()
    with open(fasta, 'rU') as fh_in:
        header = None
        seq = ""
        for line in fh_in:
            line = line.rstrip()
            if line:
                if line[0] == '>':
                    if header:
                        header2seq[header] = seq
                    header = line[1:]
                    seq = ""
                else:
                    seq += line

        if header:
            header2seq[header] = seq

    return header2seq

# # Only handles the start positions for now
# def convert_1nuc_to_1codon(nuc_startpos_1based):
#     return nuc_startpos_1based/NUC_PER_CODON + 1
#
# # Only handles the start positions for now
# def convert_1codon_to_1nuc(codon_startpos_1based):
#     return ((codon_startpos_1based - 1) * NUC_PER_CODON) + 1

def get_tree_len_depth(treefilename):
    """
    Returns tuple of (sum of all branch lengths in tree (excluding root branch), deepest root to tip distance)
    :param treefilename:
    :return: total branch length sum of tree (excluding root branch), deepest root to tip distance
    :rtype: (float, float)
    """
    tree = Phylo.read(treefilename, "newick")

    root_branch_length = tree.clade.branch_length  # set to 1.0  for some odd reason
    tree_branch_length  = tree.clade.total_branch_length()  # sum of all branch lengths including root branch length
    unroot_tree_len = tree_branch_length  - root_branch_length
    clade_depths = tree.depths()
    longest_depth = 0.0

    for clade, depth in clade_depths.iteritems():
        if clade.is_terminal():
            if longest_depth < depth:
                longest_depth = depth

    return unroot_tree_len, longest_depth