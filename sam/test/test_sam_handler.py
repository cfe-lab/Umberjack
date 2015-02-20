import unittest
import sam.sam_handler

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)

    def test_create_msa_slice_from_sam(self):
        mapping_cutoff = 20
        read_qual_cutoff = 20
        max_prop_N = 0.1
        breadth_thresh = 0.875
        start_pos = 2041
        end_pos = 2340
        is_insert = False
        ref_len = None
        ref = None

        sam_filename = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/aln/small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.consensus.bwa.sort.query.sam"
        out_fasta_filename = "./test_small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.consensus.bwa.{}_{}.fasta".format(start_pos, end_pos)
        written = sam.sam_handler.create_msa_slice_from_sam(sam_filename, ref, out_fasta_filename, mapping_cutoff, read_qual_cutoff,
                                  max_prop_N, breadth_thresh, start_pos, end_pos, is_insert, ref_len)
        print (written)

        sam_filename = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/reads/small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.sort.query.sam"
        out_fasta_filename = "./test_small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.{}_{}.fasta".format(start_pos, end_pos)
        written = sam.sam_handler.create_msa_slice_from_sam(sam_filename, ref, out_fasta_filename, mapping_cutoff, read_qual_cutoff,
                                  max_prop_N, breadth_thresh, start_pos, end_pos, is_insert, ref_len)
        print (written)
if __name__ == '__main__':
    unittest.main()
