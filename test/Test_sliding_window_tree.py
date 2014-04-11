import unittest
import sliding_window_tree

# For now, use a *.remap.sam file (paired end reads aligned to a consensus sequence with indels removed).
SAM_FILENAME = "./data/TestSample-RT_S17.HIV1B-vif.remap.sam"
MAPQ_CUTOFF = 0  # alignment quality cutoff
MAX_PROP_N = 0  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./data/TestSample-RT_S17.HIV1B-vif.10.conseq"

MIN_WINDOW_BREADTH_COV_FRACTION = 0.5
MIN_WINDOW_DEPTH_COV = 2

class TestSlidingWindowTree(unittest.TestCase):

    def test_process_windows(self):
        sliding_window_tree.process_windows(sam_filename=SAM_FILENAME,
                                            ref_fasta_filename=REFERENCE_FASTA,
                                            mapping_cutoff=MAPQ_CUTOFF,
                                            read_qual_cutoff=READ_QUAL_CUTOFF,
                                            max_prop_N=MAX_PROP_N,
                                            window_breadth_thresh=MIN_WINDOW_BREADTH_COV_FRACTION,
                                            window_depth_thresh=MIN_WINDOW_DEPTH_COV)


if __name__ == '__main__':
    unittest.main()
