import unittest
import Utility
import os


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
        expected_total_codons_by_pos = [1, 2, 0, 0, 0]
        total_codons_by_pos = Utility.get_total_codons_by_pos(msa_fasta_filename=msa_fasta_filename)

        self.assertEqual(0, cmp(total_codons_by_pos, expected_total_codons_by_pos),
                         msg="Expected total codons by positions=" + ",".join(str(x) for x in expected_total_codons_by_pos) +
                        " but got=" + ",".join(str(x) for x in total_codons_by_pos) )

        os.remove(msa_fasta_filename)


if __name__ == '__main__':
    unittest.main()
