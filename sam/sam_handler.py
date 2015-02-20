import re
import os
import sys
import subprocess
import logging
import math
from sam_constants import SamFlag as SamFlag
from sam_constants import  SamHeader as SamHeader
import sam_record
import paired_records
import Utility



NEWICK_NAME_RE = re.compile('[:;\-\(\)\[\]]')


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
console_handler.setFormatter(formatter)
LOGGER.addHandler(console_handler)









# TODO:  make stop codon remove optional
def __write_seq(fh_out, name, seq, max_prop_N, breadth_thresh):
    """
    Helper function to write out sequence to fasta file handle if has sufficient bases.
    :return  True if sequence written out
    :rtype : boolean
    :param fh_out:
    :param name:
    :param seq:
    """
    # Set stop codons to NNN so that HyPhy doesn't auto-remove a site without telling us
    for nuc_pos in range(0, len(seq), Utility.NUC_PER_CODON):
        codon = seq[nuc_pos:nuc_pos+Utility.NUC_PER_CODON]
        if Utility.CODON2AA.get(codon, "") == Utility.STOP_AA:
            seq = seq[0:nuc_pos] + "NNN" + seq[nuc_pos+Utility.NUC_PER_CODON:]

    if seq.count('N') / float(len(seq)) <= max_prop_N and (seq.count("N") + seq.count("-"))/float(len(seq)) <= (1.0-breadth_thresh):
        # Newick tree formats don't like special characters.  Convert them to underscores.
        newick_nice_qname = re.sub(pattern=NEWICK_NAME_RE, repl='_', string=name)
        fh_out.write(">" + newick_nice_qname + "\n")
        fh_out.write(seq + "\n")
        return True
    return False


# TODO:  make stop codon removal optional
def create_msa_slice_from_sam(sam_filename, ref, out_fasta_filename, mapping_cutoff, read_qual_cutoff,
                              max_prop_N, breadth_thresh, start_pos=None, end_pos=None, is_insert=False, ref_len=None):
    """
    Parse SAM file contents for query-ref (pairwise) aligned sequences for a specific reference contig.
    For paired-end reads, merges the reads into a single sequence with gaps with respect to the reference.
    Creates a pseudo-multiple sequence alignment on all the query sequences and reference.
    Writes the MSA sequences to out_fasta_filename.
    Converts query names so that they are compatible with Newick format in phylogenetic reconstruction by
        converting colons, semicolons, parentheses to underscores.


    Sam file must be query sorted.  Secondary alignments are ignored.
    NB: Only takes the primary alignment.

    :param str sam_filename: full path to sam file.  Must be queryname sorted.
    :param str ref: name of reference contig to form MSA alignments to.  If None, then splits out all alignments to any reference.
    :param str out_fasta_filename: full path of fasta file to write to.  Will completely overwite file.
    :param int mapping_cutoff:  Ignore alignments with mapping quality lower than the cutoff.
    :param int read_qual_cutoff: Convert bases with quality lower than this cutoff to N unless both mates agree.
    :param float max_prop_N:  Do not output merged sequences with proportion of N higher than the cutoff
    :param int start_pos: 1-based start nucleotide start position of slice.  If None, then uses beginning of ref.
    :param int end_pos: 1-based end nucleotide start position of slice.  If None, then uses end of ref.
    :param int ref_len: length of reference.  If None, then takes length from sam headers.
    """


    LOGGER.debug("About to slice fasta " + out_fasta_filename + " from " + sam_filename)
    if os.path.exists(out_fasta_filename) and os.path.getsize(out_fasta_filename):
        LOGGER.warn("Found existing Sliced MSA-Fasta " + out_fasta_filename + ". Not regenerating.")
        total_seq = Utility.get_total_seq_from_fasta(out_fasta_filename)
        LOGGER.debug("Done slice fasta " + out_fasta_filename)
        return total_seq

    if not is_query_sort(sam_filename):
        raise ValueError("Sam file must be queryname sorted and header must specify sort order")

    total_written = 0
    if not ref_len:
        ref_len = get_reflen(sam_filename, ref)
    with open(sam_filename, 'r') as sam_fh, open(out_fasta_filename, 'w') as out_fasta_fh:
        prev_mate = None
        is_expect_mate = False
        for line in sam_fh:
            if line.startswith(SamHeader.TAG_HEADER_PREFIX):  # skip the headers
                continue

            lines_arr = line.rstrip().split('\t')
            if len(lines_arr) < 11:  # in case there are no alignments
                break

            qname, flag_str, rname, pos_str, mapq_str, cigar, rnext, pnext, tlen, seq, qual = lines_arr[:11]
            flag = int(flag_str)
            pos = int(pos_str)
            mapq = int(mapq_str)


            # If read failed to map or has poor mapping quality, skip it
            # Be careful!  The rname is only set to '*' and pos is only set to '0'
            #   if all reads in the mate pair are unmapped.
            # If the mate is mapped, but this read is not, the rname and pos will be set to the mate's reference hit.
            # From SAM specs:
            # "For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read
            #   should have RNAME and POS identical to its mate"
            if (SamFlag.IS_UNMAPPED & flag or SamFlag.IS_SECONDARY_ALIGNMENT & flag or SamFlag.IS_CHIMERIC_ALIGNMENT & flag or
                    rname == '*' or cigar == '*' or pos == 0 or mapq < mapping_cutoff or (ref and rname != ref)):
                continue

            mate = sam_record.SamRecord(ref_len=ref_len)
            mate.fill_record_vals(qname, flag, rname, seq, cigar, mapq, qual, pos, rnext, pnext)

            is_mate_paired = SamFlag.IS_PAIRED & flag and not SamFlag.IS_MATE_UNMAPPED & flag and (rnext == "=" or not ref or rnext == ref)

            # We are expecting this record or future records to be the next mate in pair
            if is_expect_mate:
                if (not prev_mate or not SamFlag.IS_PAIRED & prev_mate.flag or SamFlag.IS_MATE_UNMAPPED & prev_mate.flag or
                         not (prev_mate.rnext == "=" or not ref or prev_mate.rnext == ref)):
                    if prev_mate:
                        LOGGER.error("Invalid logic.  Prevmate.qname=" + prev_mate.qname + " flag=" + str(prev_mate.flag) + " rnext=" + prev_mate.rnext)
                    raise ValueError("Invalid logic.  Shouldn't get here. - " + line)

                # Check if last sam record and this sam record are for the same paired read
                if prev_mate.qname == qname and is_mate_paired:
                    pair = paired_records.PairedRecord(prev_mate, mate)
                    # TODO:  get a merge_sam_reads functioin that doesn't do all the stats
                    mseq, stats = pair.merge_sam_reads(q_cutoff=read_qual_cutoff,
                                                       pad_space_btn_mates="N",
                                                       do_insert_wrt_ref=is_insert,
                                                       do_pad_wrt_ref=False,
                                                       do_pad_wrt_slice=True,
                                                       slice_start_wrt_ref_1based=start_pos,
                                                       slice_end_wrt_ref_1based=end_pos)
                    is_written = __write_seq(out_fasta_fh, pair.get_name(), mseq, max_prop_N, breadth_thresh)
                    total_written += 1 if is_written else 0
                else:
                    # TODO:  dont spit out this warning if we are missing a mate because we rejected its map qual
                    LOGGER.warn("Sam record inconsistent.  Expected pair for " + prev_mate.qname + " but got " + qname)
                    # This sam record does not pair with the previous sam record.
                    # Do low quality masking on the previous record and write it out.
                    mseq, mqual, stats = prev_mate.get_seq_qual(do_pad_wrt_ref=False,
                                                         do_pad_wrt_slice=True,
                                                         do_mask_low_qual=True,
                                                         q_cutoff=read_qual_cutoff,
                                                         slice_start_wrt_ref_1based=start_pos,
                                                         slice_end_wrt_ref_1based=end_pos,
                                                         do_insert_wrt_ref=is_insert)
                    is_written = __write_seq(out_fasta_fh, prev_mate.qname, mseq, max_prop_N, breadth_thresh)
                    total_written += 1 if is_written else 0

            if not is_mate_paired:
                # Do low quality masking on the this record and write it out.
                mseq, mqual, stats = mate.get_seq_qual(do_pad_wrt_ref=False,
                                                do_pad_wrt_slice=True,
                                                do_mask_low_qual=True,
                                                q_cutoff=read_qual_cutoff,
                                                slice_start_wrt_ref_1based=start_pos,
                                                slice_end_wrt_ref_1based=end_pos,
                                                do_insert_wrt_ref=is_insert)
                is_written = __write_seq(out_fasta_fh, mate.qname, mseq, max_prop_N, breadth_thresh)
                total_written += 1 if is_written else 0
                is_expect_mate = False
                prev_mate = None
            else:
                if is_expect_mate and prev_mate.qname == qname and is_mate_paired:  # already wrote the pair
                    is_expect_mate = False
                    prev_mate = None
                else:
                    is_expect_mate = True
                    prev_mate = mate



        # In case we went through the full sam file but couldn't find the next paired mate to map well to the same ref
        # Do low quality masking on the previous record and write it out.
        if is_expect_mate and prev_mate:
            mseq, mqual = prev_mate.get_seq_qual(do_pad_wrt_ref=False,
                                                 do_pad_wrt_slice=True,
                                                 do_mask_low_qual=True,
                                                 q_cutoff=read_qual_cutoff,
                                                 slice_start_wrt_ref_1based=start_pos,
                                                 slice_end_wrt_ref_1based=end_pos,
                                                 do_insert_wrt_ref=is_insert)
            is_written = __write_seq(out_fasta_fh, prev_mate.qname, mseq, max_prop_N, breadth_thresh)
            total_written += 1 if is_written else 0

    LOGGER.debug("Done slice fasta " + out_fasta_filename)
    return total_written



def is_query_sort(sam_filename):
    """
    :param sam_filename:
    :return: whether the sam is sorted by queryname using the header.  False if header if absent.
    :rtype: boolean
    """
    with open(sam_filename, 'rU') as fh_in:
        header = fh_in.next()  # header should be the first line if present
        header = header.rstrip()
        if not header.startswith(SamHeader.TAG_HEADER_START):
            return False
        tags = header.split(SamHeader.TAG_SEP)
        for tag in tags[1:]:
            key, val = tag.split(SamHeader.TAG_KEY_VAL_SEP)
            if key == SamHeader.TAG_SORT_ORDER_KEY:
                if val == SamHeader.TAG_SORT_ORDER_VAL_QUERYNAME:
                    return True
                else:
                    return False
    return False

def get_reflen(sam_filename, ref):
    with open(sam_filename, 'rU') as fh_in:
        for line in fh_in:
            if not line.startswith(SamHeader.TAG_HEADER_PREFIX):
                return None
            if line.startswith(SamHeader.TAG_REFSEQ_START):
                line = line.rstrip()
                is_found_ref = False
                length = None
                for tag in line.split(SamHeader.TAG_SEP)[1:]:
                    key, val = tag.split(SamHeader.TAG_KEY_VAL_SEP)
                    if key == SamHeader.TAG_REFSEQ_NAME_KEY and val == ref:
                        is_found_ref = True
                    if key == SamHeader.TAG_REFSEQ_LEN_KEY:
                        length = int(val)
                if is_found_ref:
                    return length
    return None



def create_depth_file_from_bam(bam_filename):
    """
    Gets the coverage from samtools depth.
    Creates a samtools per-base depth file with the same name as bam_filename but appended with ".depth".

    TODO: what to do with STDERR

    :param str bam_filename:  full path to sorted and indexed bam file
    :return: returns full filepath to depth file
    :rtype : str
    :raise subprocess.CalledProcessError
    """

    # Get per-base depth
    depth_filename = bam_filename + ".depth"
    with open(depth_filename, 'w') as depth_fh:
        subprocess.check_call(['samtools', 'depth', bam_filename], stdout=depth_fh, shell=False)
    return depth_filename


def sam_to_sort_bam(sam_filename, ref_filename):
    """
    Creates index of reference.  This creates a <ref_filename>.fai file.
    Converts sam to bam file sorted by coordinates.
    This creates a <sam_filename prefix>.bam and <sam filename prefix>.bam.sort files.
    Creates index of sorted bam.  This creates a <bam_filename>.index file.
    Pipes STDERR to STDOUT.
    Uses default samtools from PATH environment variable.

    :return: full filepath to sorted, indexed bam  file
    :rtype : str
    :param str sam_filename: full filepath to sam alignment file
    :param str ref_filename:  full filepath to reference fasta
    """

    sam_filename_prefix, file_suffix = os.path.splitext(sam_filename)
    bam_filename = sam_filename_prefix + ".bam"

    # index reference
    index_ref_filename = ref_filename + ".fai"
    subprocess.check_call(['samtools', 'faidx', ref_filename], stderr=subprocess.STDOUT, shell=False)

    # convert sam to bam
    subprocess.check_call(['samtools', 'view', '-Sb', '-o', bam_filename, sam_filename],
                          stderr=subprocess.STDOUT, shell=False)

    # Sort the bam file by leftmost position on the reference assembly.  Required for samtools depth.
    sorted_bam_filename_prefix = sam_filename_prefix + ".sort"
    sorted_bam_filename = sorted_bam_filename_prefix + ".bam"
    subprocess.check_call(['samtools', 'sort', bam_filename, sorted_bam_filename_prefix], stderr=subprocess.STDOUT, shell=False)

    # index sorted bam file
    subprocess.check_call(['samtools', 'index', sorted_bam_filename, ], stderr=subprocess.STDOUT, shell=False)

    return sorted_bam_filename


def get_ave_coverage_from_bam (bam_filename, ref, pos_start, pos_end):
    """
    If the depth file does not exist, then create one.
    Reads the depth file and returns average coverage in the specified region.

    :param str bam_filename: full path to bam file of alignments against the ref
    :param str ref: name of reference sequence to get coverage for
    :param int pos_start: region starting position, 1-based
    :param int pos_end: region ending position, 1-based
    :return: average per-base coverage for the specified reference and region.
    :rtype : float
    """
    import os

    depth_filename = bam_filename + ".depth"
    if not os.path.isfile(depth_filename) or os.path.getsize(depth_filename) <= 0:
        create_depth_file_from_bam(bam_filename)

    total_coverage = 0
    with open(depth_filename, 'r') as depth_fh:
        for line in depth_fh:
            # Output of samtools depth:
            # <reference>   <position>  <depth>
            line_ref, line_pos, line_depth = line.rstrip().split('\t')
            if line_ref == ref and int(line_pos) >= pos_start and int(line_pos) <= pos_end:
                total_coverage += int(line_depth)

    ave_coverage = round(float(total_coverage) / (pos_end - pos_start + 1), 0)
    return ave_coverage


def get_ref_len(sam_filename, ref):
    with open(sam_filename, 'rU') as fh:
        for line in fh:
            if line.startswith("@SQ"):
                fields = line.split("\t")
                if line.find("SN:" + ref) >= 0:
                    line_reflen_match = re.search(r"\tLN:(\d+)", line)
                    return int(line_reflen_match.group(1))

    return None