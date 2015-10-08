import unittest
import Utility
import os
import tempfile
import math


class MyTestCase(unittest.TestCase):


    def test_get_sitelist_unambig_aa(self):
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
        total_codons_by_pos = Utility.get_sitelist_unambig_aa(msa_fasta_filename=tmp_msa_fasta.name)

        self.assertEqual(0, cmp(total_codons_by_pos, expected_total_codons_by_pos),
                         msg="Expected total codons by positions=" + ",".join(str(x) for x in expected_total_codons_by_pos) +
                        " but got=" + ",".join(str(x) for x in total_codons_by_pos) )

        os.remove(tmp_msa_fasta.name)


    def test_get_depth_incomplete_codon(self):
        # The last codon is incomplete
        tmp_msa_fasta = tempfile.NamedTemporaryFile(delete=False)
        #ACT AC- NNN NCG T-T A
        #XXX GG- AC- TCT --A C
        #XN- GG- ACN --N --A A
        tmp_msa_fasta.write(">test1\n")
        tmp_msa_fasta.write("ACTAC-NNNNCGT-TA\n")
        tmp_msa_fasta.write(">test2\n")
        tmp_msa_fasta.write("---GG-AC-TCT--AC\n")
        tmp_msa_fasta.write(">test3\n")
        tmp_msa_fasta.write("-N-GG-ACN--N--AA\n")
        tmp_msa_fasta.flush()
        os.fsync(tmp_msa_fasta.file.fileno())
        tmp_msa_fasta.close()

        aln = Utility.Consensus()
        aln.parse(msa_fasta_filename=tmp_msa_fasta.name)
        expected_aln_width = 16
        actual_aln_width = aln.get_alignment_len()
        self.assertEqual(expected_aln_width, aln.get_alignment_len(),
                         "Expected alignment length=" + str(expected_aln_width) + " but got " + str(actual_aln_width))


        # Only unambiguous nucleotides allowed
        expected_nuc_depth_unambig =   [1,1,1,  3,3,0,  2,2,0,  1,2,2,  1,0,3,  3]
        expected_codon_depth_unambig = [1,      0,      0,      1,      0,      0]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
            self.assertEqual(expected_nuc_depth_unambig[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_unambig[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                self.assertEqual(expected_codon_depth_unambig[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_unambig[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))

        # Allow mixtures
        #ACT AC- NNN NCG T-T A
        #XXX GG- AC- TCT --A C
        #XN- GG- ACN --N --A A
        expected_nuc_depth_ambig =   [1,2,1,  3,3,0,  3,3,2,  2,2,3,  1,0,3,  3]
        expected_codon_depth_ambig = [1,      0,      2,      2,      0,      0]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=True, is_count_gaps=False, is_count_pad=False)
            self.assertEqual(expected_nuc_depth_ambig[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_ambig[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=True, is_count_gaps=False, is_count_pad=False)
                self.assertEqual(expected_codon_depth_ambig[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_ambig[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))


        # Allow internal gaps
        #ACT AC- NNN NCG T-T A
        #XXX GG- AC- TCT --A C
        #XN- GG- ACN --N --A A
        expected_nuc_depth_gap =   [1,1,2,  3,3,3,  2,2,1,  2,3,2,  3,3,3,  3]
        expected_codon_depth_gap = [1,      3,      1,      1,      3,      0]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=False, is_count_gaps=True, is_count_pad=False)
            self.assertEqual(expected_nuc_depth_gap[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_gap[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=False, is_count_gaps=True, is_count_pad=False)
                self.assertEqual(expected_codon_depth_gap[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_gap[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))


        # Allow external gaps
        #ACT AC- NNN NCG T-T A
        #XXX GG- AC- TCT --A C
        #XN- GG- ACN --N --A A
        expected_nuc_depth_pad =   [3,2,2,  3,3,0,  2,2,0,  1,2,2,  1,0,3,  3]
        expected_codon_depth_pad = [2,      0,      0,      1,      0,      3]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=True)
            self.assertEqual(expected_nuc_depth_pad[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_pad[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=True)
                self.assertEqual(expected_codon_depth_pad[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_pad[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))

        os.remove(tmp_msa_fasta.name)


    def test_get_depth_complete_codon(self):
        # The last codon is incomplete
        tmp_msa_fasta = tempfile.NamedTemporaryFile(delete=False)
        #ACT AC- NNN NCG T-T ACT
        #XXX GG- AC- TCT --A XXX
        #XN- GG- ACN --N --A -NX
        tmp_msa_fasta.write(">test1\n")
        tmp_msa_fasta.write("ACTAC-NNNNCGT-TACT\n")
        tmp_msa_fasta.write(">test2\n")
        tmp_msa_fasta.write("---GG-AC-TCT--A---\n")
        tmp_msa_fasta.write(">test3\n")
        tmp_msa_fasta.write("-N-GG-ACN--N--A-N-\n")
        tmp_msa_fasta.flush()
        os.fsync(tmp_msa_fasta.file.fileno())
        tmp_msa_fasta.close()

        aln = Utility.Consensus()
        aln.parse(msa_fasta_filename=tmp_msa_fasta.name)
        expected_aln_width = 18
        actual_aln_width = aln.get_alignment_len()
        self.assertEqual(expected_aln_width, aln.get_alignment_len(),
                         "Expected alignment length=" + str(expected_aln_width) + " but got " + str(actual_aln_width))


        # Only unambiguous nucleotides allowed
        expected_nuc_depth_unambig =   [1,1,1,  3,3,0,  2,2,0,  1,2,2,  1,0,3,  1,1,1]
        expected_codon_depth_unambig = [1,      0,      0,      1,      0,      1]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
            self.assertEqual(expected_nuc_depth_unambig[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_unambig[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                self.assertEqual(expected_codon_depth_unambig[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_unambig[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))

        # Allow mixtures
        #ACT AC- NNN NCG T-T ACT
        #XXX GG- AC- TCT --A XXX
        #XN- GG- ACN --N --A -NX
        expected_nuc_depth_ambig =   [1,2,1,  3,3,0,  3,3,2,  2,2,3,  1,0,3,  1,2,1]
        expected_codon_depth_ambig = [1,      0,      2,      2,      0,      1]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=True, is_count_gaps=False, is_count_pad=False)
            self.assertEqual(expected_nuc_depth_ambig[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_ambig[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=True, is_count_gaps=False, is_count_pad=False)
                self.assertEqual(expected_codon_depth_ambig[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_ambig[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))


        # Allow internal gaps
        #ACT AC- NNN NCG T-T ACT
        #XXX GG- AC- TCT --A XXX
        #XN- GG- ACN --N --A -NX
        expected_nuc_depth_gap =   [1,1,2,  3,3,3,  2,2,1,  2,3,2,  3,3,3,  2,1,1]
        expected_codon_depth_gap = [1,      3,      1,      1,      3,      1]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=False, is_count_gaps=True, is_count_pad=False)
            self.assertEqual(expected_nuc_depth_gap[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_gap[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=False, is_count_gaps=True, is_count_pad=False)
                self.assertEqual(expected_codon_depth_gap[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_gap[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))


        # Allow external gaps
        #ACT AC- NNN NCG T-T ACT
        #XXX GG- AC- TCT --A XXX
        #XN- GG- ACN --N --A -NX
        expected_nuc_depth_pad =   [3,2,2,  3,3,0,  2,2,0,  1,2,2,  1,0,3,  2,2,3]
        expected_codon_depth_pad = [2,      0,      0,      1,      0,      2]
        for nucpos in range(expected_aln_width):
            actual_nuc_depth = aln.get_depth(pos_0based=nucpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=True)
            self.assertEqual(expected_nuc_depth_pad[nucpos], actual_nuc_depth,
                             "Expected nuc depth=" + str(expected_nuc_depth_pad[nucpos]) + " but got " +
                             str(actual_nuc_depth) + " at nucpos=" + str(nucpos))

            if nucpos % 3 == 0:
                codonpos = nucpos/3
                actual_codon_depth = aln.get_codon_depth(codon_pos_0based=codonpos, is_count_ambig=False, is_count_gaps=False, is_count_pad=True)
                self.assertEqual(expected_codon_depth_pad[codonpos], actual_codon_depth,
                             "Expected codon depth=" + str(expected_codon_depth_pad[codonpos]) + " but got " +
                             str(actual_codon_depth) + " at codonpos=" + str(codonpos))

        os.remove(tmp_msa_fasta.name)



    def test_get_consensus_from_msa(self):
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.write(">seq1\n")
        tmpfile.write("ATNG\n")
        tmpfile.write(">seq2\n")
        tmpfile.write("GCGC\n")
        tmpfile.write(">seq3\n")
        tmpfile.write("GTGC\n")
        tmpfile.flush()
        os.fsync(tmpfile.file.fileno())
        tmpfile.close()
        cons = Utility.Consensus()
        cons.parse(tmpfile.name)
        cons_str = cons.get_consensus()
        self.assertEqual("GTGC", cons_str)


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

        expected = [ -(1.0/3 * math.log(1.0/3, 2) + 2.0/3 * math.log(2.0/3, 2))/math.log(2, 2),
                     -(2.0/3 * math.log(2.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(2, 2),
                     0,
                     -(1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(3, 2)]

        for pos in range (0, len(expected)):
            actual = cons.get_metric_entropy(pos)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)


        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.write(">seq1\n")
        tmpfile.write("ATNG\n")
        tmpfile.write(">seq2\n")
        tmpfile.write("GCGT\n")
        tmpfile.write(">seq3\n")
        tmpfile.write("GTGC\n")
        tmpfile.write(">seq4\n")
        tmpfile.write("NC--\n")
        tmpfile.flush()
        os.fsync(tmpfile.file.fileno())
        tmpfile.close()

        cons = Utility.Consensus()
        cons.parse(tmpfile.name)

        expected = [ -(1.0/3 * math.log(1.0/3, 2) + 2.0/3 * math.log(2.0/3, 2))/math.log(2, 2),
                     -(2.0/4 * math.log(2.0/4, 2) + 2.0/4 * math.log(2.0/4, 2))/math.log(2, 2),
                     0,
                     -(1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(3, 2)]

        for pos in range (0, len(expected)):
            actual = cons.get_metric_entropy(pos)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)


    def test_metric_entropy_ambig(self):
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

        expected = [ -(1.0/3 * math.log(1.0/3, 2) + 2.0/3 * math.log(2.0/3, 2))/math.log(2, 2),
                     -(2.0/3 * math.log(2.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(2, 2),
                     -(2.25/3 * math.log(2.25/3, 2) + 3*0.25/3 * math.log(0.25/3, 2))/math.log(4, 2),
                     -(1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(3, 2)]

        for pos in range (0, len(expected)):
            actual = cons.get_metric_entropy(pos, is_count_ambig=True)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)


    def test_metric_entropy_gap(self):
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.write(">seq1\n")
        tmpfile.write("AT-G\n")
        tmpfile.write(">seq2\n")
        tmpfile.write("GCGT\n")
        tmpfile.write(">seq3\n")
        tmpfile.write("GTGC\n")
        tmpfile.flush()
        os.fsync(tmpfile.file.fileno())
        tmpfile.close()

        cons = Utility.Consensus()
        cons.parse(tmpfile.name)

        expected = [ -(1.0/3 * math.log(1.0/3, 2) + 2.0/3 * math.log(2.0/3, 2))/math.log(2, 2),
                     -(2.0/3 * math.log(2.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(2, 2),
                     -(2.25/3 * math.log(2.25/3, 2) + 3*0.25/3 * math.log(0.25/3, 2))/math.log(4, 2),
                     -(1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2) + 1.0/3 * math.log(1.0/3, 2))/math.log(3, 2)]

        for pos in range (0, len(expected)):
            actual = cons.get_metric_entropy(pos, is_count_ambig=True, is_count_gaps=True)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)


    def test_get_conserve(self):
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

        expected = [ 2/3.0, 2/3.0, 2.0/2, 1.0/3]

        for pos in range (0, len(expected)):
            actual = cons.get_conserve(pos)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)

    def test_get_conserve_ambig(self):
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

        expected = [ 2/3.0, 2/3.0, 2.25/3, 1.0/3]

        for pos in range (0, len(expected)):
            actual = cons.get_conserve(pos, is_count_ambig=True)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)

    def test_get_conserve_gap(self):
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.write(">seq1\n")
        tmpfile.write("AT-G\n")
        tmpfile.write(">seq2\n")
        tmpfile.write("GCGT\n")
        tmpfile.write(">seq3\n")
        tmpfile.write("GTGC\n")
        tmpfile.flush()
        os.fsync(tmpfile.file.fileno())
        tmpfile.close()

        cons = Utility.Consensus()
        cons.parse(tmpfile.name)

        expected = [ 2/3.0, 2/3.0, 2.25/3, 1.0/3]

        for pos in range (0, len(expected)):
            actual = cons.get_conserve(pos, is_count_ambig=True, is_count_gaps=True)
            print("Pos={} expected={} actual={}".format(pos, expected[pos], actual))
            self.assertEqual(expected[pos], actual, "Pos={} expected={} actual={}".format(pos, expected[pos], actual))

        os.remove(tmpfile.name)


if __name__ == '__main__':
    unittest.main()
