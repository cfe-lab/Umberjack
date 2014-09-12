import os
import Utility
import csv
import glob

# Columns in the HyPhy dN/dS tab-separates values file
HYPHY_TSV_DN_COL = 'dN'
HYPHY_TSV_DS_COL = 'dS'
HYPHY_TSV_S_COL = 'Observed S Changes'
HYPHY_TSV_N_COL = 'Observed NS Changes'
HYPHY_TSV_SCALED_DN_MINUS_DS_COL = 'Scaled dN-dS'
HYPHY_TSV_PROB_FULLSEQ_S_COL = 'P{S leq. observed}'
HYPHY_TSV_NEG_PROB_FULLSEQ_S_COL = 'P{S geq. observed}'
HYPHY_TSV_EXP_S_COL = 'E[S Sites]'
HYPHY_TSV_EXP_N_COL = 'E[NS Sites]'

class _SiteDnDsInfo:
    """
    Internal use only.  Keeps track of selection information at a codon site.
    """
    def __init__(self):
        self.accum_win_dnds = 0.0
        self.accum_win_dn_minus_ds = 0.0
        self.total_win_cover_site = 0
        self.total_syn_subs = 0
        self.total_nonsyn_subs = 0
        self.total_reads = 0
        self.total_exp_syn_subs = 0
        self.total_exp_nonsyn_subs = 0
        self.window_dnds_subs = [] # list of number of substitutions

    def add_dnds(self, dnds, dn_minus_ds, reads, syn_subs, nonsyn_subs, exp_syn_subs, exp_nonsyn_subs):
        """
        Insert selection information from a window for the codon site.
        :param dnds: dN/dS for this codon site from a single window.
        :param dn_minus_ds: scaled dN-dS for this codon site from a single window.
        :param reads: total reads that contain a valid codon (not - or N's) at this codon site in the window.
        :param syn_subs: total synonymoous substitutions for this codon site in the window
        :param nonsyn_subs: total nonsynonymous substitutions for this codon site in the window
        :param exp_syn_subs:  expected number of synonymous substitutions
        :param exp_nonsyn_subs: expected nonsynomous substitutions
        """
        self.total_win_cover_site += 1
        if dnds is not None:
            self.accum_win_dnds += dnds
        self.accum_win_dn_minus_ds += dn_minus_ds
        self.total_reads += reads
        self.total_syn_subs += syn_subs
        self.total_nonsyn_subs += nonsyn_subs
        self.total_exp_syn_subs += exp_syn_subs
        self.total_exp_nonsyn_subs += exp_nonsyn_subs
        self.window_dnds_subs.append([dnds, syn_subs + nonsyn_subs, reads])

    def get_weighted_ave_dnds(self):
        """
        Return weighted average dN/dS from all windows for the codon site.
        Average weighted by number of substitutions for the codon site in the window.
        :rtype : float
        """
        if not self.window_dnds_subs:
            return None
        else:
            total_weighted_site_dnds_over_all_windows = 0.0
            total_subs_at_site_over_all_windows = self.total_syn_subs + self.total_nonsyn_subs
            for (site_dnds, subs_at_site_in_window, reads) in self.window_dnds_subs:
                total_weighted_site_dnds_over_all_windows += (site_dnds * subs_at_site_in_window)
            return total_weighted_site_dnds_over_all_windows / total_subs_at_site_over_all_windows

    def get_weighted_byreads_ave_dnds(self):
        """
        Return weighted average dN/dS from all windows for the codon site.
        Average weighted by number of substitutions for the codon site in the window.
        :rtype : float
        """
        if not self.window_dnds_subs:
            return None
        else:
            total_weighted_site_dnds_over_all_windows = 0.0
            for (site_dnds, subs_at_site_in_window, reads_in_window) in self.window_dnds_subs:
                total_weighted_site_dnds_over_all_windows += (site_dnds * reads_in_window)
            return total_weighted_site_dnds_over_all_windows / self.total_reads


    def get_ave_dnds(self):
        """
        Return average dN/dS from all windows for the codon site
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return self.accum_win_dnds / self.total_win_cover_site

    def get_ave_dn_minus_ds(self):
        """
        Return average scaled dN-dS from all windows for the codon site
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return self.accum_win_dn_minus_ds / self.total_win_cover_site

    def get_window_coverage(self):
        """
        Return total windows that cover this codon site
        :rtype : int
        """
        return self.total_win_cover_site

    def get_ave_read_coverage(self):
        """
        Return average reads that cover this codon site with a valid codon (is not - or Ns) over all windows.
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return float(self.total_reads)/self.total_win_cover_site

    def get_ave_syn_subs(self):
        """
        Return average synonymous substitutions at this codon site over all windows.
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return float(self.total_syn_subs)/self.total_win_cover_site

    def get_ave_nonsyn_subs(self):
        """
        Return average nonsynonymous substitutions at this codon site over all windows.
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return float(self.total_nonsyn_subs)/self.total_win_cover_site

    def get_ave_subs(self):
        """
        Return average substitutions at this codon site over all windows.
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return float(self.total_nonsyn_subs + self.total_syn_subs)/self.total_win_cover_site


class SeqDnDsInfo:
    """
    Keeps track of selection information at a sequence of codon sites.
    """
    def __init__(self, seq_len):
        """
        Constructor.
        :rtype : SeqDnDsInfo
        :param int seq_len: total codons in the sequence
        """
        self.dnds_seq = [_SiteDnDsInfo() for i in range(seq_len)]

    def add_site_dnds(self, site_1based, dnds, dn_minus_ds, reads, syn_subs, nonsyn_subs, exp_syn_subs, exp_nonsyn_subs):
        """
        Keep track of dn/ds from a window containing this site.
        :param site_1based: 1-based codon site
        :param dnds: dn/ds for this codon site as determined from a window fed into hyphy
        :param dn_minus_ds:
        :param reads:
        :param syn_subs:
        :param nonsyn_subs:
        :param exp_syn_subs:  expected synonymous substitutions
        :param exp_nonsyn_subs: expected nonsynomous substitutions
        """
        self.dnds_seq[site_1based-1].add_dnds(dnds=dnds, dn_minus_ds=dn_minus_ds, reads=reads,
                                              syn_subs=syn_subs, nonsyn_subs=nonsyn_subs,
                                              exp_syn_subs=exp_syn_subs, exp_nonsyn_subs=exp_nonsyn_subs)

    def get_site_ave_dnds(self, site_1based):
        """
        Get the average dN/dS at the codon site over all windows.

        :return: average dN/dS at the codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_ave_dnds()

    def get_site_ave_dn_minus_ds(self, site_1based):
        """
        Get the average scaled dN - dS at the codon site over all windows.

        :return: average scaled dN - dS at the codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_ave_dn_minus_ds()

    def get_seq_len(self):
        """
        Get the length of the codon sequence.

        :return: length of codon sequence
        :rtype : int
        """
        return len(self.dnds_seq)

    def get_site_window_cov(self, site_1based):
        """
        Get the total windows that cover the given codon site.

        :return: total windows that cover the given codon site.
        :rtype : int
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_window_coverage()

    def get_site_ave_read_cov(self, site_1based):
        """
        Get the average number of reads that cover the given codon site over all windows.
        The reads must contain a valid codon at the given codon site (as opposed to - or N)

        :return: average number of reads that cover the given codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_ave_read_coverage()

    def get_site_ave_syn_subs(self, site_1based):
        """
        Get the average number of synonymous substitutions at the given codon site over all windows.

        :return: average number of synonymous substitutions at the given codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_ave_syn_subs()

    def get_site_ave_nonsyn_subs(self, site_1based):
        """
        Get the average number of nonsynonymous substitutions at the given codon site over all windows.

        :return: average number of nonsynonymous substitutions at the given codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_ave_nonsyn_subs()

    def get_site_ave_subs(self, site_1based):
        """
        Get the average number of substitutions at the given codon site over all windows.

        :return: average number of substitutions at the given codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_ave_subs()

    def get_weighted_site_ave_dnds(self, site_1based):
        """
        Get the average number of substitutions at the given codon site over all windows,
        weighted by the number of substitutions at the site in a window

        :return: average number of substitutions at the given codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_weighted_ave_dnds()

    def get_weighted_byreads_ave_dnds(self, site_1based):
        """
        Get the average number of substitutions at the given codon site over all windows,
        weighted by the number of substitutions at the site in a window

        :return: average number of substitutions at the given codon site over all windows.
        :rtype : float
        :param int site_1based : 1-based codon site position
        """
        return self.dnds_seq[site_1based-1].get_weighted_byreads_ave_dnds()

    def get_multisite_weighted_ave_dnds(self, site_start_1based, site_end_1based):
        """
        Gets the average dN/dS weighted by substiutions across a range of codon sites over all windows covering the codon range
        :param site_start_1based:
        :param site_end_1based:
        :return:
        """
        total_dnds = 0.0
        if site_start_1based > site_end_1based:
            raise ValueError("Codon range start must be before range end")

        total_sites = site_end_1based - site_start_1based + 1
        for site_0based in range(site_start_1based-1, site_end_1based):
            total_dnds += self.dnds_seq[site_0based].get_weighted_ave_dnds()

        return total_dnds / total_sites

    def get_multisite_ave_dnds(self, site_start_1based, site_end_1based):
        """
        Gets the average dN/dS across a range of codon sites over all windows covering the codon range.
        :param site_start_1based:
        :param site_end_1based:
        :return:
        """
        total_dnds = 0.0
        if site_start_1based > site_end_1based:
            raise ValueError("Codon range start must be before range end")

        total_nonsyn = 0.0
        total_syn = 0.0
        total_exp_nonsyn = 0.0
        total_exp_syn = 0.0


        for site_0based in range(site_start_1based-1, site_end_1based):
            total_nonsyn += self.dnds_seq[site_0based].total_nonsyn_subs
            total_syn += self.dnds_seq[site_0based].total_syn_subs
            total_exp_nonsyn += self.dnds_seq[site_0based].total_exp_nonsyn_subs
            total_exp_syn += self.dnds_seq[site_0based].total_exp_syn_subs

        if total_syn == 0 and total_nonsyn == 0:
            return None
        else:
            return (total_nonsyn / total_exp_nonsyn) * (total_exp_syn / total_syn)


def __get_window_seq(seq, start_pos, end_pos, breadth_thresh=0.0):
    """
    If the sequence fits the window constraints, then returns the sliced sequence.
    Otherwise returns None.

    :return:  Sequence that fits in the window, right padded if necessary.  None if sequence doesn't fit window constraints.
    :rtype : str
    :param str seq: sequence with no newlines
    :param int start_pos : 1-based start position of window
    :param int end_pos : 1-based end position of window
    :param float breadth_thresh : fraction of sequence that be A,C,T,or G within start_pos and end_pos inclusive.
    """
    if seq and len(seq) >= start_pos:
        slice_seq = seq[start_pos-1:end_pos]
        slice_seq_uc = slice_seq.upper()
        a_count = slice_seq_uc.count('A')
        t_count = slice_seq_uc.count('T')
        g_count = slice_seq_uc.count('G')
        c_count = slice_seq_uc.count('C')
        if float(a_count+g_count+c_count+t_count)/(end_pos - start_pos + 1) >= breadth_thresh:
            pad = '-' * max(0, end_pos - len(seq))
            return slice_seq

    return None


def create_slice_msa_fasta(fasta_filename, out_fasta_filename, start_pos, end_pos, breadth_thresh=0.0):
    """
    From a fasta file of multiple sequence alignments, extract the sequence sequence from desired region.
    Writes the extracted sequences into a new fasta file with ".<start_pos>_<end_pos>.fasta" suffix.
    If the sequence is shorter than <end_pos>, then fills in the gaps with '-' characters so that it ends at <end_pos>

    Only puts in the read into the sliced msa fasta if it obeys the window constraints.

    :return:  total sequences written
    :rtype : int
    :param str fasta_filename: full file path to fasta of multiple sequence alignments
    :param str out_fasta_filename:  full file path to output fasta of sliced multiple sequence alignments
    :param int start_pos : 1-based start position of region to extract
    :param int end_pos: 1-based end position of region to extract
    :param float breadth_thresh: fraction of sequence that be A,C,T,or G within start_pos and end_pos inclusive.
    """

    total_seq = 0
    with open(fasta_filename, 'r') as fasta_fh, open(out_fasta_filename, 'w') as slice_fasta_fh:
        header = ""
        seq = ""
        for line in fasta_fh:
            line = line.rstrip()
            if line:
                line = line.split()[0]  # remove trailing whitespace and any test after the first whitespace

                if line[0] == '>':  # previous sequence is finished.  Write out previous sequence
                    window_seq = __get_window_seq(seq=seq, start_pos=start_pos, end_pos=end_pos, breadth_thresh=breadth_thresh)
                    if window_seq:
                        slice_fasta_fh.write(header + "\n")
                        slice_fasta_fh.write(window_seq + "\n")
                        total_seq += 1

                    seq = ""
                    header = line

                else:   # cache current sequence so that entire sequence is on one line
                    seq += line

        window_seq = __get_window_seq(seq=seq, start_pos=start_pos, end_pos=end_pos, breadth_thresh=breadth_thresh)
        if window_seq:
            slice_fasta_fh.write(header + "\n")
            slice_fasta_fh.write(window_seq + "\n")
            total_seq += 1

    return total_seq


def get_seq_dnds_info(dnds_tsv_dir, pvalue_thresh, ref, ref_codon_len):
    """
    Get dN/dS information from multiple windows for multiple codon sites as a SeqDnDsInfo object.
    We expect that all HyPhy has written out dn/ds tsv files for every window.

    :return: a SeqDnDsInfo object containing the information about selection at all codon sites from multiple windows
    :rtype: SeqDnDsInfo
    :param str dnds_tsv_dir:  directory containing sitewise dn/ds tab-separated files generated by HyPhy
    :param float pvalue_thresh:  pvalue threshold for selection significantly different from neutral evolution
    :param str ref: reference contig/chromosome name.  Does not include the '>' faster header delimtier or
                    any thing after the first space in the header text
    :param int ref_codon_len:  length of the reference in codons

    Assumes these are the columns in the HyPhy DN/DS TSV output:
    - Observed S Changes
    - Observed NS Changes
    - E[S Sites]: proportion of random one-nucleotide substitutions that are expected to be synonymous
    - E[NS Sites]: proportion of random one-nucleotide substitutions that are expected to be non-synonymous
    - Observed S. Prop.: observed proportion of synomymous substitutions = Observed S Changes / (Observed S Changes + Observed NS Changes)
    - P{S}:  proportion of substitutions expected to be synonymous under neutral evolution = E[S Sites]/(E[S Sites] + E[NS Sites])
    - dS: observed synonymous substitutions / expected proportion synonymous substitutions = Observed S Changes / E[S Sites]
    - dN: observed non synonymous substitutions / expected proportion nonsynonymous substitutions = Observed NS Changes / E[NS Sites]
    - dN-dS:  difference between dS and dN
    - P{S leq. observed}:  binomial distro pvalue.  Probability of getting less than the observed synonymous substitutions
        under the binomial distribution where probability of 1 synonymous codon = P{S}
    - P{S geq. observed}: binomial distro  pvalue (for the other tail).  Probability of getting more than the observed synonymous substitutions
        under the binomial distribution where probability of 1 synonymous codon = P{S}
    - Scaled dN-dS:  dN-dS normalized by the total length of the tree.

    We only want significant Dn/Ds.
    """
    seq_dnds_info = SeqDnDsInfo(seq_len=ref_codon_len)

    for dnds_tsv_filename in glob.glob(dnds_tsv_dir + os.sep + "*.dnds.tsv"):
        with open(dnds_tsv_filename, 'r') as dnds_fh:
            # *.{start bp}_{end bp}.dnds.tsv filenames use 1-based nucleotide position numbering
            dnds_tsv_fileprefix = dnds_tsv_filename.split('.dnds.tsv')[0]
            win_nuc_range = dnds_tsv_fileprefix.split('.')[-1]
            # Window ends at this 1-based nucleotide position with respect to the reference
            win_start_nuc_pos_1based_wrt_ref = int(win_nuc_range.split('_')[0])
            # Window starts at this 1-based codon position with respect to the reference
            win_start_codon_1based_wrt_ref = win_start_nuc_pos_1based_wrt_ref/Utility.NUC_PER_CODON + 1

            msa_slice_fasta_filename = dnds_tsv_fileprefix + ".fasta"
            codons_by_window_pos = Utility.get_total_codons_by_pos(msa_fasta_filename=msa_slice_fasta_filename)

            reader = csv.DictReader(dnds_fh, delimiter='\t',)
            for offset, codon_row in enumerate(reader):    # Every codon site is a row in the *.dnds.tsv file
                pval = float(codon_row[HYPHY_TSV_PROB_FULLSEQ_S_COL])
                neg_tail_pval = float(codon_row[HYPHY_TSV_NEG_PROB_FULLSEQ_S_COL])
                if pval <= pvalue_thresh or neg_tail_pval <= pvalue_thresh:
                    dN = float(codon_row[HYPHY_TSV_DN_COL])
                    dS = float(codon_row[HYPHY_TSV_DS_COL])
                    dn_minus_ds = float(codon_row[HYPHY_TSV_SCALED_DN_MINUS_DS_COL])
                    syn_subs = float(codon_row[HYPHY_TSV_S_COL])
                    nonsyn_subs = float(codon_row[HYPHY_TSV_N_COL])
                    exp_syn_subs = float(codon_row[HYPHY_TSV_EXP_S_COL])
                    exp_nonsyn_subs = float(codon_row[HYPHY_TSV_EXP_N_COL])
                    ref_codon_1based = win_start_codon_1based_wrt_ref + offset
                    codons = codons_by_window_pos[offset]

                    if dS == 0:
                        dnds = None
                    else:
                        dnds = dN/dS
                    seq_dnds_info.add_site_dnds(site_1based=ref_codon_1based, dnds=dnds, dn_minus_ds=dn_minus_ds,
                                                reads=codons, syn_subs=syn_subs, nonsyn_subs=nonsyn_subs,
                                                exp_syn_subs=exp_syn_subs, exp_nonsyn_subs=exp_nonsyn_subs)

    return seq_dnds_info


def tabulate_dnds(dnds_tsv_dir, ref, ref_nuc_len, pvalue_thresh, output_dnds_tsv_filename, comments):
    """
    Aggregate selection information from multiple windows for each codon site.
    Output selection information into a tab separated file with the following columns:
    - Ref:  reference name
    - Site:  1-based codon site in the reference
    - dNdS:  dN/dS
    - dN_minus_dS:  dN - dS  (this has not been scaled by the number of substitutions at the codon site)
    - Codons: the number of reads containing a valid codon at this codon site
    - NonSyn: the number of nonsynonymous substitutions at this codon site amongst the reads
    - Syn:  the number of synonymous substitutions at this codon site amongst the reads
    - Subst:  the number of substitutions at this codon site amongst the reads

    :return: a SeqDnDsInfo object containing the information about selection at all codon sites from all windows
    :rtype : SeqDnDsInfo
    :param str dnds_tsv_dir:  output directory of dN/dS tab separated files generated by HyPhy
    :param str ref: name of reference contig
    :param int ref_nuc_len:  length of reference contig in nucleotides
    :param float pvalue_thresh: p-value threshold for significant selection
    :param str output_dnds_tsv_filename: full filepath of aggregated selection tsv to write to
    :param str comments: any comments to add at the top of the aggregated selection tsv
    """
    seq_dnds_info = get_seq_dnds_info(dnds_tsv_dir=dnds_tsv_dir, pvalue_thresh=pvalue_thresh, ref=ref,
                                                    ref_codon_len=ref_nuc_len/Utility.NUC_PER_CODON)

    SMOOTH_DIST = 15
    with open(output_dnds_tsv_filename, 'w') as dnds_fh:
        dnds_fh.write("# " + comments + "\n")
        dnds_fh.write("Ref\tSite\tdNdSWeightBySubst\tdN_minus_dS\tWindows\tCodons\tNonSyn\tSyn\tSubst\tdNdSWeightByReads\tmultisitedNdS\n")
        for site in range(1, seq_dnds_info.get_seq_len() + 1):
            site_dnds = seq_dnds_info.get_weighted_site_ave_dnds(site_1based=site)
            site_dn_minus_ds = seq_dnds_info.get_site_ave_dn_minus_ds(site_1based=site)
            window = seq_dnds_info.get_site_window_cov(site_1based=site)
            reads = seq_dnds_info.get_site_ave_read_cov(site_1based=site)
            nonsyn = seq_dnds_info.get_site_ave_nonsyn_subs(site_1based=site)
            syn = seq_dnds_info.get_site_ave_syn_subs(site_1based=site)
            subs = seq_dnds_info.get_site_ave_subs(site_1based=site)
            site_dnds_weight_by_reads = seq_dnds_info.get_weighted_byreads_ave_dnds(site_1based=site)

            smooth_dist_start = max(site-SMOOTH_DIST, 1)
            smooth_dist_end = min(site+SMOOTH_DIST, seq_dnds_info.get_seq_len())
            multisite_dnds = seq_dnds_info.get_multisite_ave_dnds(site_start_1based=smooth_dist_start, site_end_1based=smooth_dist_end)

            line = "\t".join([ref, str(site), str(site_dnds), str(site_dn_minus_ds), str(window), str(reads),  str(nonsyn), str(syn), str(subs),
                              str(site_dnds_weight_by_reads),
                              str(multisite_dnds)])
            dnds_fh.write(line + "\n")
    return seq_dnds_info