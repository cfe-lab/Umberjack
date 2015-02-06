import unittest
import Utility
import os
import tempfile
import math

MSA_FASTA_FILENAME = "./simulations/data/sample_genomes.fas"
CONSENSUS_FASTA_FILENAME = "./simulations/out/test_consensus.fasta"

class MyTestCase(unittest.TestCase):


    def test_get_total_pos_seq_from_fasta(self):
        tmp_msa_fasta = tempfile.NamedTemporaryFile(delete=False)
        tmp_msa_fasta.write(">test1\n")
        #ACT ACN NNN NCG T-T
        tmp_msa_fasta.write("ACTACNNNNNCGT-TA\n")
        tmp_msa_fasta.write(">test2\n")
        # --- GG- -G- --N --A
        tmp_msa_fasta.write("---GG--G---N--AC\n")
        tmp_msa_fasta.write(">test3\n")
        # --- GG- -G- --N --A AA
        tmp_msa_fasta.write("---GG--G---N--AAAT\n")
        tmp_msa_fasta.flush()
        os.fsync(tmp_msa_fasta.file.fileno())
        tmp_msa_fasta.close()

        expected_total_codons_by_pos = [1, 3, 0, 0, 0, 1]
        total_codons_by_pos = Utility.get_total_codons_by_pos(msa_fasta_filename=tmp_msa_fasta.name)

        self.assertEqual(0, cmp(total_codons_by_pos, expected_total_codons_by_pos),
                         msg="Expected total codons by positions=" + ",".join(str(x) for x in expected_total_codons_by_pos) +
                        " but got=" + ",".join(str(x) for x in total_codons_by_pos) )

        os.remove(tmp_msa_fasta.name)

    def test_get_consensus_from_msa(self):
        # TODO:  make trivial example instead of this
        pass


    def test_shannon_entropy(self):
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.write(">seq1\n")
        tmpfile.write("ATNG\n")
        tmpfile.write(">seq2\n")
        tmpfile.write("GCGT\n")
        tmpfile.write(">seq3\n")
        tmpfile.write("GTGC\n")
        tmpfile.flush()
        os.fsync(tmpfile.file.fileno())
        tmpfile.close()

        cons = Utility.Consensus()
        cons.parse(tmpfile.name)

        expected = [ -(1.0/3 * math.log(1.0/3, 2) + 2.0/3 * math.log(2.0/3, 2)),
                     -(2.0/3 * math.log(2.0/3, 2) + 1.0/3 * math.log(1.0/3, 2)),
                     -(2.0/2 * math.log(2.0/2, 2) ),
                     -(1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))]

        for pos in range (0, len(expected)):
            actual = cons.get_shannon_entropy(pos)
            self.assertEqual(expected[pos], actual, "Pos=0 expected={} actual={}".format(expected[pos], actual))

        os.remove(tmpfile.name)


    def test_metric_entropy(self):
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.write(">seq1\n")
        tmpfile.write("ATNG\n")
        tmpfile.write(">seq2\n")
        tmpfile.write("GCGT\n")
        tmpfile.write(">seq3\n")
        tmpfile.write("GTGC\n")
        tmpfile.flush()
        os.fsync(tmpfile.file.fileno())
        tmpfile.close()

        cons = Utility.Consensus()
        cons.parse(tmpfile.name)

        expected = [ -(1.0/3 * math.log(1.0/3, 2) + 2.0/3 * math.log(2.0/3, 2))/3,
                     -(2.0/3 * math.log(2.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/3,
                     -(2.0/2 * math.log(2.0/2, 2) )/3,
                     -(1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/3]

        for pos in range (0, len(expected)):
            actual = cons.get_metric_entropy(pos)
            print("Pos=0 expected={} actual={}".format(expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos=0 expected={} actual={}".format(expected[pos], actual))

        os.remove(tmpfile.name)

if __name__ == '__main__':
    unittest.main()
