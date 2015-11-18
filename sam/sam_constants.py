import re

PHRED_SANGER_OFFSET = 33
CIGAR_RE = re.compile('[0-9]+[MIDNSHPX=]')
SEQ_PAD_CHAR = '-'
QUAL_PAD_CHAR = ' '     # This is the ASCII character right below lowest PHRED quality score in Sanger qualities  (-1)
QUAL_ZERO_CHAR = str(chr(PHRED_SANGER_OFFSET + 0))  # Quality of zero
SAM_UNSPECIFIED = "*"

class SamFlag:
    IS_PAIRED =                0x001
    IS_MAPPED_IN_PROPER_PAIR = 0x002
    IS_UNMAPPED =              0x004
    IS_MATE_UNMAPPED =         0x008
    IS_REVERSE =               0x010
    IS_MATE_REVERSE =          0x020
    IS_FIRST =                 0x040
    IS_SECOND =                0x080
    IS_SECONDARY_ALIGNMENT =   0x100
    IS_FAILED =                0x200
    IS_DUPLICATE =             0x400
    IS_CHIMERIC_ALIGNMENT =    0x800

class SamHeader:
    TAG_HEADER_PREFIX = "@"
    TAG_HEADER_START = "@HD"
    TAG_SEP = "\t"
    TAG_KEY_VAL_SEP = ":"
    TAG_SORT_ORDER_KEY = "SO"
    TAG_SORT_ORDER_VAL_QUERYNAME = "queryname"
    TAG_REFSEQ_START = "@SQ"
    TAG_REFSEQ_NAME_KEY = "SN"
    TAG_REFSEQ_LEN_KEY = "LN"

# Mandatory fields in a sam file
SAM_FIELDS = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
SAM_FIELD_DELIM = "\t"

