import unittest
import Utility
import os


MSA_FASTA_FILENAME = "./simulations/data/sample_genomes.fas"
CONSENSUS_FASTA_FILENAME = "./simulations/data/out/test_consensus.fasta"

class MyTestCase(unittest.TestCase):


    def test_get_total_pos_seq_from_fasta(self):
        msa_fasta_filename = "./data/out/utility_fake_msa.fasta"
        with open(msa_fasta_filename, 'w') as fh:
            fh.write(">test1\n")
            #ACT ACN NNN NCG T-T
            fh.write("ACTACNNNNNCGT-T\n")
            fh.write(">test2\n")
            # --- GG- -G- --N --A
            fh.write("---GG--G---N--A\n")
            fh.write(">test3\n")
            # --- GG- -G- --N --A AA
            fh.write("---GG--G---N--AAA\n")
        expected_total_codons_by_pos = [1, 3, 0, 0, 0, 1]
        total_codons_by_pos = Utility.get_total_codons_by_pos(msa_fasta_filename=msa_fasta_filename)

        self.assertEqual(0, cmp(total_codons_by_pos, expected_total_codons_by_pos),
                         msg="Expected total codons by positions=" + ",".join(str(x) for x in expected_total_codons_by_pos) +
                        " but got=" + ",".join(str(x) for x in total_codons_by_pos) )

        os.remove(msa_fasta_filename)

    def test_get_consensus_from_msa(self):
        Utility.get_consensus_from_msa(MSA_FASTA_FILENAME, CONSENSUS_FASTA_FILENAME)

        # TODO: do a diff on CONSENSUS_FASTA_FILENAME and /home/thuy/gitrepo/SlidingWindow/test/simulations/data/sample_genomes.consensus.fasta


if __name__ == '__main__':
    unittest.main()
