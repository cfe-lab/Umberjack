import Utility
import csv
import glob
import os
import hyphy.hyphy_handler as hyphy_handler
import Bio.SeqIO as SeqIO
import re
import math

def collect_dnds(output_dir, output_csv_filename, full_popn_fasta, comments=None):
    """
    Collects everything related to dnds into 1 table.  Does not do any aggregation of values.  Useful for debugging.
    :return:
    """
    with open(output_csv_filename, 'w') as fh_out:
        if comments:
            fh_out.write(comments)

        writer = csv.DictWriter(fh_out, fieldnames=["Window_Start", "Window_End",
                                                    "Reads",  # Max read depth for the window (not necessary for the codon site)
                                                    "CodonSite",  # 1-based codon site
                                                    "CodonDepth",  # Total unambiguous codon (depth) at the codon site
                                                    "Conserve",  # Average per-base fraction of conservation across the codon
                                                    "Entropy",  # Average per-base metric entropy across the codon
                                                    "N",  # Observed Nonsynonymous substitutions
                                                    "S",  # Observed Nonsynonymous substitutions
                                                    "EN",  # Expected Nonsynonymous substitutions
                                                    "ES",  # Expected Synonymous substitutions
                                                    "dN", "dS",
                                                    "dN_minus_dS", # dN-dS scaled by the total substitutions at the codon site
                                                    "AmbigBase",  # N nucleotide
                                                    "Pad", # left or right pad gap
                                                    "Err",  # Nucleotide errors within the codon
                                                    "Err_N", # nonsynonymous AA change due to sequence error
                                                    "Err_S",  # synonymous AA change due to sequence error
                                                    "Ambig_N",  # Ambiguous base changes the AA.  Should be always 0
                                                    "Ambig_S", # ambigous base does not change the AA.
                                                    "TreeLen",  # Tree length
                                                    "TreeDepth", # deepest tip to root distance
                                                    ]
                                )
        writer.writeheader()
        for slice_fasta_filename in glob.glob(output_dir + os.sep + "*.fasta"):

            # *.{start bp}_{end bp}.fasta filenames use 1-based nucleotide position numbering
            slice_fasta_fileprefix = slice_fasta_filename.split('.fasta')[0]

            tree_filename = slice_fasta_filename.replace(".fasta", ".tree")
            tree_len, tree_depth = Utility.get_tree_len_depth(tree_filename)

            win_nuc_range = slice_fasta_fileprefix.split('.')[-1]
            # Window ends at this 1-based nucleotide position with respect to the reference
            win_start_nuc_pos_1based_wrt_ref, win_end_nuc_pos_1based_wrt_ref = [int(x) for x in win_nuc_range.split('_')]
            # Window starts at this 1-based codon position with respect to the reference
            win_start_codon_1based_wrt_ref = win_start_nuc_pos_1based_wrt_ref/Utility.NUC_PER_CODON + 1


            codons_by_window_pos = Utility.get_total_codons_by_pos(msa_fasta_filename=slice_fasta_filename)
            reads = Utility.get_total_seq_from_fasta(slice_fasta_filename)
            consensus = Utility.Consensus()
            consensus.parse(slice_fasta_filename)

            ns, pad, seq_err, err_aa_change, err_aa_nochange, ambig_aa_change, ambig_aa_nochange = error_by_codonpos(slice_fasta_filename, win_start_nuc_pos_1based_wrt_ref, full_popn_fasta)

            outrow = dict()
            outrow["Window_Start"] = win_start_nuc_pos_1based_wrt_ref
            outrow["Window_End"] = win_end_nuc_pos_1based_wrt_ref
            outrow["Reads"] = reads

            dnds_tsv_filename = slice_fasta_filename.replace(".fasta", ".dnds.tsv")
            if os.path.exists(dnds_tsv_filename):
                with open(dnds_tsv_filename, 'rU') as fh_dnds_tsv:
                    reader = csv.DictReader(fh_dnds_tsv, delimiter='\t')
                    for offset, codon_row in enumerate(reader):    # Every codon site is a row in the *.dnds.tsv file
                        outrow["CodonSite"] = win_start_codon_1based_wrt_ref + offset
                        outrow["CodonDepth"] = codons_by_window_pos[offset]
                        outrow["Conserve"] = consensus.get_ave_conserve(offset, offset + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                        outrow["Entropy"] = consensus.get_ave_shannon_entropy(offset, offset + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                        outrow["N"] = float(codon_row[hyphy_handler.HYPHY_TSV_N_COL])
                        outrow["S"] = float(codon_row[hyphy_handler.HYPHY_TSV_S_COL])
                        outrow["ES"] = float(codon_row[hyphy_handler.HYPHY_TSV_EXP_S_COL])
                        outrow["EN"] = float(codon_row[hyphy_handler.HYPHY_TSV_EXP_N_COL])
                        outrow["dN"] = float(codon_row[hyphy_handler.HYPHY_TSV_DN_COL])
                        outrow["dS"] = float(codon_row[hyphy_handler.HYPHY_TSV_DS_COL])
                        outrow["dN_minus_dS"] = float(codon_row[hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])

                        outrow["AmbigBase"] = ns[offset]
                        outrow["Pad"] = pad[offset]
                        outrow["Err"] = seq_err[offset]
                        outrow["Err_N"] = err_aa_change[offset]
                        outrow["Err_S"] = err_aa_nochange[offset]
                        outrow["Ambig_N"] = ambig_aa_change[offset]
                        outrow["Ambig_S"] = ambig_aa_nochange[offset]
                        outrow["TreeLen"] = tree_len
                        outrow["TreeDepth"] = tree_depth
                        writer.writerow(outrow)
            else:
                for offset, codons in enumerate(codons_by_window_pos):
                    outrow["CodonSite"] = win_start_codon_1based_wrt_ref + offset
                    outrow["CodonDepth"] = codons
                    outrow["Conserve"] = consensus.get_ave_conserve(offset, offset + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                    outrow["Entropy"] = consensus.get_ave_shannon_entropy(offset, offset + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                    outrow["AmbigBase"] = ns[offset]
                    outrow["Pad"] = pad[offset]
                    outrow["Err"] = seq_err[offset]
                    outrow["Err_N"] = err_aa_change[offset]
                    outrow["Err_S"] = err_aa_nochange[offset]
                    outrow["Ambig_N"] = ambig_aa_change[offset]
                    outrow["Ambig_S"] = ambig_aa_nochange[offset]
                    outrow["TreeLen"] = tree_len
                    outrow["TreeDepth"] = tree_depth
                    writer.writerow(outrow)


def error_by_codonpos(slice_msa_fasta, slice_start_wrt_ref_1based, full_popn_fasta):
    # Collect read error stats per window - codon

    #full_popn_recdict = SeqIO.to_dict(SeqIO.parse("/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/small.cov2.indiv1k.codon400.bwa.rand.mixed.fasta", "fasta"))
    full_popn_recdict = SeqIO.to_dict(SeqIO.parse(full_popn_fasta, "fasta"))


    longest_seq = Utility.get_longest_seq_size_from_fasta(slice_msa_fasta)
    total_codons = int(math.ceil(float(longest_seq)/Utility.NUC_PER_CODON))
    ns = [0] * total_codons
    pad = [0] * total_codons
    seq_err = [0] * total_codons
    err_aa_change = [0] * total_codons
    err_aa_nochange = [0] * total_codons
    ambig_aa_change = [0] * total_codons  # aa change due to ambiguous base.  Should be 0
    ambig_aa_nochange = [0] * total_codons  # no aa change even though ambiguous base


    for record in SeqIO.parse(slice_msa_fasta, "fasta"):
        template, read = record.id.split("_")

        # The read is a subset of the slice
        # first nongap char 0-based nucleotide position wrt slice
        read_start_wrt_slice_0based = re.search(r"[^\-]", str(record.seq)).start()
        # last nongap char 0-based nucleotide position wrt slice
        read_end_wrt_slice_0based = re.search(r"[^\-][\-]*$", str(record.seq)).start()

        # The 0-based codon position with respect to the slice for the first codon that isn't left padded with gap chars
        read_codon_start_wrt_slice_0based = int(math.ceil(float(read_start_wrt_slice_0based)/Utility.NUC_PER_CODON))
        # The 0-based codon position with respect to the slice for the last codon that isn't right padded with gap chars
        read_codon_end_wrt_slice_0based = ((read_end_wrt_slice_0based+1)/Utility.NUC_PER_CODON) - 1

        # NB:  the nonpadded portion of the read might not start on a codon, but the MSA should always start on a codon
        for nuc_pos_wrt_slice_0based in range(0, longest_seq, Utility.NUC_PER_CODON):
            read_codon = str(record.seq[nuc_pos_wrt_slice_0based:nuc_pos_wrt_slice_0based + Utility.NUC_PER_CODON])
            read_codon += "-" * (Utility.NUC_PER_CODON - len(read_codon))  # right-pad with gaps in case the nucleotide sequence doesn't end on codon end

            nuc_pos_wrt_ref_0based = slice_start_wrt_ref_1based-1 + nuc_pos_wrt_slice_0based
            template_codon = str(full_popn_recdict[template].seq[nuc_pos_wrt_ref_0based:nuc_pos_wrt_ref_0based + Utility.NUC_PER_CODON])

            codon_pos_wrt_slice_0based = nuc_pos_wrt_slice_0based/Utility.NUC_PER_CODON
            is_codon_has_err = False
            is_codon_has_ambig = False
            for i in range(0, Utility.NUC_PER_CODON):
                if read_codon[i] == "N":
                    is_codon_has_ambig = True
                    ns[codon_pos_wrt_slice_0based] += 1
                elif (read_codon[i] == "-" and (codon_pos_wrt_slice_0based < read_codon_start_wrt_slice_0based or
                                                        codon_pos_wrt_slice_0based > read_codon_end_wrt_slice_0based)):
                    pad[codon_pos_wrt_slice_0based] += 1
                    is_codon_has_ambig = True
                elif read_codon[i] != template_codon[i]:  # sequence error
                    is_codon_has_err = True
                    seq_err[codon_pos_wrt_slice_0based] += 1

            # aa change or no aa-change due to sequence error
            if is_codon_has_err and Utility.CODON2AA.get(read_codon) and Utility.CODON2AA.get(template_codon):
                if Utility.CODON2AA.get(read_codon) == Utility.CODON2AA.get(template_codon):
                    err_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    err_aa_change[codon_pos_wrt_slice_0based] += 1
            elif is_codon_has_ambig and Utility.CODON2AA.get(read_codon) and Utility.CODON2AA.get(template_codon):
                if Utility.CODON2AA.get(read_codon) == Utility.CODON2AA.get(template_codon):
                    ambig_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    ambig_aa_change[codon_pos_wrt_slice_0based] += 1


    return ns, pad, seq_err, err_aa_change, err_aa_nochange, ambig_aa_change, ambig_aa_nochange
