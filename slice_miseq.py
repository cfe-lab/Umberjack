import os
import Utility
import csv
import glob
import logging
import fasttree.fasttree_handler as fasttree
import re
import hyphy.hyphy_handler as hyphy_handler

LOGGER = logging.getLogger(__name__)



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
        self.total_reads_nolowsyn = 0.0
        self.total_dn_minus_ds_nolowsyn = 0.0
        self.total_dn_minus_ds_weightby_reads_nolowsyn = 0.0
        self.total_dnds_nolowsyn = 0.0
        self.accum_win_dnds_weightby_reads_nolowsyn = 0.0

    def add_dnds(self, dnds, dn_minus_ds, reads, syn_subs, nonsyn_subs, exp_syn_subs, exp_nonsyn_subs):
        """
        Insert selection information from a window for the codon site.
        :param dnds: dN/dS for this codon site from a single window.
        :param dn_minus_ds: scaled dN-dS for this codon site from a single window.
        :param reads: total reads that contain a valid codon (not - or N's) at this codon site in the window.
        :param syn_subs: total synonymous substitutions for this codon site in the window
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
        if syn_subs >= 1.0:
            self.total_reads_nolowsyn += reads
            self.total_dn_minus_ds_nolowsyn += dn_minus_ds
            self.total_dn_minus_ds_weightby_reads_nolowsyn += (reads*dn_minus_ds)
            self.accum_win_dnds_weightby_reads_nolowsyn += (reads*dnds)

        self.window_dnds_subs.append([dnds, syn_subs, nonsyn_subs, reads, exp_syn_subs, exp_nonsyn_subs])


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
            total_subs_at_site_over_all_windows = 0
            for (site_dnds, site_syn_subst, site_nonsyn_subst, reads, exp_syn_subs, exp_nonsyn_subs) in self.window_dnds_subs:
                if not site_dnds is None:   # sometimes there are zero observed synonymous substitutions.  Thus dn/ds = None.
                    total_weighted_site_dnds_over_all_windows += (site_dnds * (site_syn_subst+site_nonsyn_subst))
                    total_subs_at_site_over_all_windows += (site_syn_subst+site_nonsyn_subst)
            if not total_subs_at_site_over_all_windows:
                return None
            return total_weighted_site_dnds_over_all_windows / total_subs_at_site_over_all_windows


    def get_ave_dnds_weightby_reads(self):
        """
        Return weighted average dN/dS from all windows for the codon site.
        Average weighted by number of non-gap reads for the codon site in the window.
        :rtype : float
        """
        if not self.window_dnds_subs:
            return None
        else:
            total_weighted_site_dnds_over_all_windows = 0.0
            total_reads_at_site_over_all_windows = 0
            for (site_dnds, site_syn_subst, site_nonsyn_subst, reads_in_window, exp_syn_subs, exp_nonsyn_subs) in self.window_dnds_subs:
                # sometimes there are zero observed synonymous substitutions, thus dn/ds = None.  Don't include those in average dn/ds
                if not site_dnds is None:
                    total_weighted_site_dnds_over_all_windows += (site_dnds * reads_in_window)
                    total_reads_at_site_over_all_windows += reads_in_window
            # Don't use self.total_reads because that counts the reads even for window-sites where dn/ds = None
            if not total_reads_at_site_over_all_windows:
                return None
            return total_weighted_site_dnds_over_all_windows / total_reads_at_site_over_all_windows


    def get_ave_dnds_weightby_reads_nolowsyn_aveall(self):
        """
        Returns weighted average dN/dS from all windows for the codon site.
        Excludes any windows in which the number of synonymous substitutions < 1.
        Also does summed of each N, EN, S, ES before calcaulting end dN/dS instead of ave of dN/dS.
        :return float:
        """
        if not self.window_dnds_subs:
            return None
        else:
            total_reads_at_site_over_all_windows = 0.0
            total_nonsyn_subs_at_site = 0.0
            total_exp_nonsyn_subs_at_site = 0.0
            total_syn_subs_at_site = 0.0
            total_exp_syn_subs_at_site = 0.0
            for (site_dnds, site_syn_subst, site_nonsyn_subst, reads_in_window, exp_syn_subs, exp_nonsyn_subs) in self.window_dnds_subs:
                # sometimes there are less than one observed synonymous substitutions.  Don't include those in average dn/ds
                # These are pretty iffy
                if site_syn_subst >= 1.0:
                    total_nonsyn_subs_at_site += (site_nonsyn_subst * reads_in_window)
                    total_exp_nonsyn_subs_at_site += (exp_nonsyn_subs * reads_in_window)
                    total_syn_subs_at_site += (site_syn_subst * reads_in_window)
                    total_exp_syn_subs_at_site += (exp_syn_subs * reads_in_window)
                    total_reads_at_site_over_all_windows += reads_in_window

            # Don't use self.total_reads because that counts the reads for every window
            if not total_reads_at_site_over_all_windows:
                return None
            # dN/dS = (total non synon subst / total exp non syn subst ) / (total syn subst / total exp syn subst)
            # = (total non synon subst * total exp syn subst ) / (total exp non syn subst * total syn subst)
            return (total_nonsyn_subs_at_site * total_exp_syn_subs_at_site)/(total_exp_nonsyn_subs_at_site * total_syn_subs_at_site )


    def get_ave_dnds_weightby_reads_nolowsyn(self):
        """
        Returns weighted average dN/dS from all windows for the codon site.
        Excludes any windows in which the number of synonymous substitutions < 1.
        :return float:
        """
        if not self.total_reads_nolowsyn:
            return None
        return self.accum_win_dnds_weightby_reads_nolowsyn/self.total_reads_nolowsyn



    def get_ave_dnds(self):
        """
        Return average dN/dS from all windows for the codon site.
        NB:  if there is dn/dS=None, then we don't want to count that window
        :rtype : float
        """
        total_site_dnds_over_all_windows = 0.0
        total_windows = 0
        # Sometimes dN/dS=None because there are 0 synonumous substitutions.  Don't include those windows
        for (site_dnds, site_syn_subst, site_nonsyn_subst, reads_in_window, exp_syn_subs, exp_nonsyn_subs) in self.window_dnds_subs:
            if not site_dnds is None:   # sometimes there are zero observed synonymous substitutions.  Thus dn/ds = None.
                    total_site_dnds_over_all_windows += site_dnds
                    total_windows += 1
        if not total_windows:
            return None
        return total_site_dnds_over_all_windows / total_windows


    def get_ave_dn_minus_ds(self):
        """
        Return average scaled dN-dS from all windows for the codon site
        :rtype : float
        """
        if not self.total_win_cover_site:
            return None
        else:
            return self.accum_win_dn_minus_ds / self.total_win_cover_site

    def get_ave_dn_minus_ds_weightby_reads_nolowsyn(self):
        """
        Return average scaled dN-dS across all windows for the codon site, weighted by reads per window.
        Ignore windows in which there are less than 1 synonymous substitution at the site
        :return float:
        """
        if not self.total_reads_nolowsyn:
            return None

        return self.total_dn_minus_ds_weightby_reads_nolowsyn/self.total_reads_nolowsyn


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
        return self.dnds_seq[site_1based-1].get_ave_dnds_weightby_reads()

    def get_ave_dnds_weightby_reads_nolowsyn(self, site_1based):
        return self.dnds_seq[site_1based-1].get_ave_dnds_weightby_reads_nolowsyn()

    def get_ave_dnds_weightby_reads_nolowsyn_aveall(self, site_1based):
        return self.dnds_seq[site_1based-1].get_ave_dnds_weightby_reads_nolowsyn_aveall()

    def get_ave_dn_minus_ds_weightby_reads_nolowsyn(self, site_1based):
        return self.dnds_seq[site_1based-1].get_ave_dn_minus_ds_weightby_reads_nolowsyn()

    def get_multisite_weighted_bysubst_ave_dnds(self, site_start_1based, site_end_1based):
        """
        Gets the average dN/dS weighted by substiutions across a range of codon sites over all windows covering the codon range
        :param site_start_1based:
        :param site_end_1based:
        :return:
        """
        total_dnds = 0.0
        if site_start_1based > site_end_1based:
            raise ValueError("Codon range start must be before range end")

        total_sites = 0.0
        for site_0based in range(site_start_1based-1, site_end_1based):
            site_weighted_ave_dnds = self.dnds_seq[site_0based].get_weighted_ave_dnds()
            if site_weighted_ave_dnds is not None:
                total_dnds += site_weighted_ave_dnds
                total_sites += 1

        if not total_sites:
            return None

        return total_dnds / total_sites

    def get_multisite_weighted_byreads_ave_dnds(self, site_start_1based, site_end_1based):
        """
        Gets the average dN/dS weighted by reads across a range of codon sites over all windows covering the codon range
        :param site_start_1based:
        :param site_end_1based:
        :return:
        """
        total_dnds = 0.0
        if site_start_1based > site_end_1based:
            raise ValueError("Codon range start must be before range end")

        total_sites = 0.0
        for site_0based in range(site_start_1based-1, site_end_1based):
            site_weighted_ave_dnds = self.dnds_seq[site_0based].get_weighted_ave_dnds()
            if site_weighted_ave_dnds is not None:
                total_dnds += site_weighted_ave_dnds
                total_sites += 1

        if not total_sites:
            return None

        return total_dnds / total_sites

    def get_multisite_ave_dnds(self, site_start_1based, site_end_1based):
        """
        Finds the average for the set of site-average dN/dS in the specified codon range.
        :param site_start_1based:
        :param site_end_1based:
        :return:
        """
        total_dnds = 0.0
        if site_start_1based > site_end_1based:
            raise ValueError("Codon range start must be before range end")

        # total_nonsyn = 0.0
        # total_syn = 0.0
        # total_exp_nonsyn = 0.0
        # total_exp_syn = 0.0
        #
        #
        # for site_0based in range(site_start_1based-1, site_end_1based):
        #     total_nonsyn += self.dnds_seq[site_0based].total_nonsyn_subs
        #     total_syn += self.dnds_seq[site_0based].total_syn_subs
        #     total_exp_nonsyn += self.dnds_seq[site_0based].total_exp_nonsyn_subs
        #     total_exp_syn += self.dnds_seq[site_0based].total_exp_syn_subs
        #
        # if total_syn == 0 or total_exp_nonsyn == 0:
        #     return None
        # else:
        #     return (total_nonsyn / total_exp_nonsyn) * (total_exp_syn / total_syn)
        total_site_ave_dnds = 0.0
        total_sites = 0
        for site_0based in range(site_start_1based-1, site_end_1based):
            site_ave_dnds = self.dnds_seq[site_0based].get_ave_dnds()
            if site_ave_dnds is not None:
                total_site_ave_dnds += site_ave_dnds
                total_sites += 1

        if total_sites == 0:
            return None
        else:
            return total_site_ave_dnds / total_sites


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



def get_seq_dnds_info(dnds_tsv_dir, ref, ref_codon_len):
    """
    Get dN/dS information from multiple windows for multiple codon sites as a SeqDnDsInfo object.
    We expect that all HyPhy has written out dn/ds tsv files for every window.

    :return: a SeqDnDsInfo object containing the information about selection at all codon sites from multiple windows
    :rtype: SeqDnDsInfo
    :param str dnds_tsv_dir:  directory containing sitewise dn/ds tab-separated files generated by HyPhy
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
                dN = float(codon_row[hyphy_handler.HYPHY_TSV_DN_COL])
                dS = float(codon_row[hyphy_handler.HYPHY_TSV_DS_COL])
                dn_minus_ds = float(codon_row[hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])
                syn_subs = float(codon_row[hyphy_handler.HYPHY_TSV_S_COL])
                nonsyn_subs = float(codon_row[hyphy_handler.HYPHY_TSV_N_COL])
                exp_syn_subs = float(codon_row[hyphy_handler.HYPHY_TSV_EXP_S_COL])
                exp_nonsyn_subs = float(codon_row[hyphy_handler.HYPHY_TSV_EXP_N_COL])
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




def tabulate_dnds(dnds_tsv_dir, ref, ref_nuc_len, output_csv_filename, comments):
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
    :param str output_csv_filename: full filepath of aggregated selection tsv to write to
    :param str comments: any comments to add at the top of the aggregated selection tsv
    """
    seq_dnds_info = get_seq_dnds_info(dnds_tsv_dir=dnds_tsv_dir, ref=ref,
                                      ref_codon_len=ref_nuc_len / Utility.NUC_PER_CODON)


    with open(output_csv_filename, 'w') as dnds_fh:
        dnds_fh.write("# " + comments + "\n")
        writer = csv.DictWriter(dnds_fh,
                                fieldnames=["Ref", "Site", "aveDnDs", "dNdSWeightBySubst", "dN_minus_dS", "Windows",
                                            "Codons", "NonSyn", "Syn", "Subst", "dNdSWeightByReads",
                                            "multisiteAvedNdS", "multisitedNdSWeightBySubst", "dNdSWeightByReadsNoLowSyn",
                                            "dNdSWeightByReadsNoLowSynAveAll", "dnMinusDsWeightByReadsNoLowSyn"])

        writer.writeheader()
        for site in range(1, seq_dnds_info.get_seq_len() + 1):
            outrow = {}
            outrow["Ref"] = ref
            outrow["Site"] = site
            outrow["aveDnDs"] = seq_dnds_info.get_site_ave_dnds(site_1based=site)
            outrow["dNdSWeightBySubst"] = seq_dnds_info.get_weighted_site_ave_dnds(site_1based=site)
            outrow["dN_minus_dS"] = seq_dnds_info.get_site_ave_dn_minus_ds(site_1based=site)
            outrow["Windows"] = seq_dnds_info.get_site_window_cov(site_1based=site)
            outrow["Codons"] = seq_dnds_info.get_site_ave_read_cov(site_1based=site)
            outrow["NonSyn"] = seq_dnds_info.get_site_ave_nonsyn_subs(site_1based=site)
            outrow["Syn"] = seq_dnds_info.get_site_ave_syn_subs(site_1based=site)
            outrow["Subst"] = seq_dnds_info.get_site_ave_subs(site_1based=site)
            outrow["dNdSWeightByReads"] = seq_dnds_info.get_weighted_byreads_ave_dnds(site_1based=site)
            outrow["dNdSWeightByReadsNoLowSyn"] = seq_dnds_info.get_ave_dnds_weightby_reads_nolowsyn(site_1based=site)
            outrow["dNdSWeightByReadsNoLowSynAveAll"] = seq_dnds_info.get_ave_dnds_weightby_reads_nolowsyn_aveall(site_1based=site)
            outrow["dnMinusDsWeightByReadsNoLowSyn"] = seq_dnds_info.get_ave_dn_minus_ds_weightby_reads_nolowsyn(site_1based=site)
            writer.writerow(outrow)

    return seq_dnds_info


def tabulate_nuc_subst(nucmodelfit_dir, output_csv_filename, comments):
    # Hyphy creates a *.nucmodelfit file that contains the best fit model (according to AIC) with this entry.  Parse it.
    #       Model averaged rates relative to AG (REV estimates):
    #           AC =   0.1902	(  0.1781)
    #           AT =   0.2058	(  0.2198)
    #           CG =   0.0573	(  0.0567)
    #           CT =   1.2453	(  1.2953)
    #           GT =   0.4195	(  0.4246)
    import fnmatch
    # .../out/RunABC/HIV1B-nef/ABC_S89.HIV1B-nef.msa.1_300.nucmodelfit
    with  open(output_csv_filename,'w') as fh_nucmodelcsv:
        fh_nucmodelcsv.write("#" + comments + "\n")
        fh_nucmodelcsv.write("ID,Ref,Window_Start,Window_End,StartBase,EndBase,Mutation,Rate\n")
        for root, dirs, filenames in os.walk(nucmodelfit_dir):
            for nucmodelfit_filename in fnmatch.filter(filenames, '*.nucmodelfit'):

                # ASSUME that multiple sequence aligned file used as input for the nucleotide model fit file is in the same folder
                # TODO:  be more general
                msa_slice_fasta_filename = nucmodelfit_filename.replace(".nucmodelfit", ".fasta")
                #nongap_by_window_pos = Utility.get_total_nongap_nuc_by_pos(msa_fasta_filename=msa_slice_fasta_filename)
                with open(os.path.join(root, nucmodelfit_filename), 'r') as fh_fit:
                    is_found_rates = False
                    # TODO:  make more general
                    sample_id, ref, msa, window, ext = os.path.basename(nucmodelfit_filename).split(".")
                    window_start, window_end = window.split("_")

                    for line in fh_fit:
                        line = line.rstrip().lstrip()
                        if not len(line):
                            continue
                        if not is_found_rates and "Model averaged rates relative to AG" in line:
                            is_found_rates = True
                            continue


                        if is_found_rates:
                            if line.find("Model averaged selection") >= 0:
                                break  # end of rates

                            match = re.findall(r'([ACGT][ACGT])\s*=\s*(\d+\.\d+)\s*\(\s*(\d+\.\d+)\s*\)', line, re.IGNORECASE)
                            if not match:
                                raise ValueError("Line should contain model rates but it doesn't: " + line)
                            subst, sym_rate, nonsym_rate = match[0]  # list of 1 tuple
                            init_base, end_base = list(subst)
                            mutation = init_base + end_base
                            fh_nucmodelcsv.write(",".join(str(x) for x in [sample_id, ref,
                                                                           window_start, window_end,
                                                                           init_base, end_base, mutation, nonsym_rate]) + "\n")


def tabulate_rates(fasttree_output_dir, output_csv_filename, comments):
    """
    Collects all the GTR model rates from all the fasttree logs in a directory and puts them into output_csv_filename.
    ASSUME that multiple sequence aligned file is in the same folder
    :param output_dir:
    :return:
    """
    import fnmatch
    # .../out/RunABC/HIV1B-nef/ABC_S89.HIV1B-nef.msa.1_300.fasttree.log
    with  open(output_csv_filename,'w') as fh_out:

        fh_out.write("#" + comments + "\n")
        #writer = csv.DictWriter(fh_out, fieldnames=["ID","Ref","Window_Start","Window_End","Window_Reads","Non_Gap_Window_Start","Mutation,Rate"])
        fh_out.write("ID,Ref,Window_Start,Window_End,Window_Reads,Non_Gap_Window_Start,Mutation,Rate\n")
        for root, dirs, filenames in os.walk(fasttree_output_dir):
            for fasttree_log in fnmatch.filter(filenames, '*.fasttree.log'):
                fullpath_fasttree_log = os.path.join(root, fasttree_log)
                AC, AG, AT, CG, CT, GT = fasttree.extract_gtr_rates(fullpath_fasttree_log)
                rates = {"AC":AC, "AG":AG, "AT":AT, "CG":CG, "CT":CT, "GT":GT}

                msa_slice_fasta_filename = fullpath_fasttree_log.replace(".fasttree.log", ".fasta")
                # sample_id.ref.msa.window_start_window_end.fasta
                name_split = os.path.basename(msa_slice_fasta_filename).split(".")
                window = name_split[-2]
                ref = name_split[-4]  # TODO:  what if reference has . in it?
                sample_id = ".".join(name_split[0:-4])
                window_start, window_end = window.split("_")
                nongap_window_start = Utility.get_total_nongap_nuc_by_pos(msa_slice_fasta_filename, 0)
                reads = Utility.get_total_seq_from_fasta(msa_slice_fasta_filename)


                for mutation, rate in rates.iteritems():
                    fh_out.write(",".join([sample_id,
                                  ref,
                                  window_start,
                                  window_end,
                                  str(reads),
                                  str(nongap_window_start),
                                  mutation,
                                  str(rate)]) + "\n")