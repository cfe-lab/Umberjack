import unittest
import slice_miseq
import os
import tempfile
import csv
import hyphy_handler
import shutil
OUT_DIR = os.path.dirname(os.path.realpath(__file__)) + os.sep + "out"


class TestSliceMiSeq(unittest.TestCase):

    def setUp(self):
        if not os.path.exists(OUT_DIR):
            shutil.rmtree(OUT_DIR)
            os.makedirs(OUT_DIR)

        self.windows = [
            # Window0:  all good
            {"syn_subs": 40.0, "nonsyn_subs":30.0, "exp_syn_subs":10.0, "exp_nonsyn_subs":10.0, "treelen":100.0, "reads":40, "startcdn":1, "endcdn":2},
             # Window1:  low nonsyn
            {"syn_subs": 40.0, "nonsyn_subs":0.4, "exp_syn_subs":10.0, "exp_nonsyn_subs":10.0, "treelen":100.0, "reads":40, "startcdn":2, "endcdn":4},
            # Window2:  all good
            {"syn_subs": 50.0, "nonsyn_subs":60.0, "exp_syn_subs":10.0, "exp_nonsyn_subs":10.0, "treelen":100.0, "reads":50, "startcdn":2, "endcdn":3},
            # Window3:  low syn
            {"syn_subs": 0.1, "nonsyn_subs":50.0, "exp_syn_subs":10.0, "exp_nonsyn_subs":10.0, "treelen":100.0, "reads":50, "startcdn":2, "endcdn":3}
        ]

        # Site-Window Coverage.  Each window repeats the same values for every site it covers
        # Site  0   1   2   3
        # Win0  *   *
        # Win1      *   *   *
        # Win2      *   *
        # Win3      *   *
        # The i'th index represents which windows cover the i'th site
        self.site_2_window = [[0], [0, 1, 2, 3], [1, 2, 3], [1]]


    def test_site_dnds_info_dnds(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows

        for window in windows:
            site_dnds_info.add_window(
                dnds=(window["nonsyn_subs"]*window["exp_syn_subs"])/(window["syn_subs"]*window["exp_nonsyn_subs"]),
                dn_minus_ds=((window["nonsyn_subs"]/window["exp_nonsyn_subs"]) - (window["syn_subs"]/window["exp_syn_subs"]))/window["treelen"],
                reads=window["reads"],
                syn_subs=window["syn_subs"], nonsyn_subs=window["nonsyn_subs"],
                exp_syn_subs=window["exp_syn_subs"], exp_nonsyn_subs=window["exp_nonsyn_subs"])

        # Include all windows, weight dn/ds by subs
        exp_total_subs = sum([window["nonsyn_subs"] + window["syn_subs"] for window in windows])
        expected_ave = 0
        for window in windows:
            expected_ave += window["nonsyn_subs"]*window["exp_syn_subs"]*(window["nonsyn_subs"]+window["syn_subs"])/(window["syn_subs"]*window["exp_nonsyn_subs"]*exp_total_subs)

        actual_ave = site_dnds_info.get_ave_dnds_weightby_subs(is_exclude_low_sub=False)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn/ds weighted by subs=" + str(expected_ave) +
                                   " but got " + str(actual_ave))

        # Include all windows, weight dn/ds by reads
        exp_total_reads = sum([window["reads"] for window in windows])
        expected_ave = 0
        for window in windows:
            expected_ave += window["nonsyn_subs"]*window["exp_syn_subs"]*window["reads"]/(window["syn_subs"]*window["exp_nonsyn_subs"]*exp_total_reads)

        actual_ave = site_dnds_info.get_ave_dnds_weightby_reads(is_exclude_low_sub=False)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn/ds weighted by reads=" + str(expected_ave) +
                                   " but got " + str(actual_ave))


        # Exclude low sub windows, weight dn-ds by subs
        exp_total_subs =sum([window["nonsyn_subs"] + window["syn_subs"] for window in [windows[0], windows[2]]])
        expected_ave = 0
        for window in [windows[0], windows[2]]:
            expected_ave += window["nonsyn_subs"]*window["exp_syn_subs"]*(window["nonsyn_subs"]+window["syn_subs"])/(window["syn_subs"]*window["exp_nonsyn_subs"]*exp_total_subs)

        actual_ave = site_dnds_info.get_ave_dnds_weightby_subs(is_exclude_low_sub=True)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn/ds weighted by subs, exclude low sub sites=" + str(expected_ave) +
                                   " but got " + str(actual_ave))

        # Exclude low sub windows, weight dn/ds by subs
        exp_total_reads = sum([window["reads"] for window in [windows[0], windows[2]]])
        expected_ave = 0
        for window in [windows[0], windows[2]]:
            expected_ave += window["nonsyn_subs"]*window["exp_syn_subs"]*window["reads"]/(window["syn_subs"]*window["exp_nonsyn_subs"]*exp_total_reads)

        actual_ave = site_dnds_info.get_ave_dnds_weightby_reads(is_exclude_low_sub=True)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn/ds weighted by reads, exclude low sub sites=" + str(expected_ave) +
                                   " but got " + str(actual_ave))


    def test_site_dnds_info_dnMinusDs(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows

        for window in windows:
            site_dnds_info.add_window(
                dnds=(window["nonsyn_subs"]*window["exp_syn_subs"])/(window["syn_subs"]*window["exp_nonsyn_subs"]),
                dn_minus_ds=((window["nonsyn_subs"]/window["exp_nonsyn_subs"]) - (window["syn_subs"]/window["exp_syn_subs"]))/window["treelen"],
                reads=window["reads"],
                syn_subs=window["syn_subs"], nonsyn_subs=window["nonsyn_subs"],
                exp_syn_subs=window["exp_syn_subs"], exp_nonsyn_subs=window["exp_nonsyn_subs"])

        # Include all windows, weight dn-ds/treelen by subs
        exp_total_subs = sum([window["nonsyn_subs"] + window["syn_subs"] for window in windows])
        expected_ave = 0
        for window in windows:
            expected_ave += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"]) - (window["syn_subs"]/window["exp_syn_subs"]))*(window["nonsyn_subs"]+window["syn_subs"])/(window["treelen"]*exp_total_subs)

        actual_ave = site_dnds_info.get_ave_dn_minus_ds_weightby_subs(is_exclude_low_sub=False)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn-ds weighted by subs=" + str(expected_ave) +
                                   " but got " + str(actual_ave))

        # Include all windows, weight dn-ds/treelen by reads
        exp_total_reads = sum([window["reads"] for window in windows])
        expected_ave = 0
        for window in windows:
            expected_ave += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"]) - (window["syn_subs"]/window["exp_syn_subs"]))*(window["reads"])/(window["treelen"]*exp_total_reads)

        actual_ave = site_dnds_info.get_ave_dn_minus_ds_weightby_reads(is_exclude_low_sub=False)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn-ds weighted by reads=" + str(expected_ave) +
                                   " but got " + str(actual_ave))


        # Exclude low sub windows, weight dn-ds/treelen by subs
        exp_total_subs =sum([window["nonsyn_subs"] + window["syn_subs"] for window in [windows[0], windows[2]]])
        expected_ave = 0
        for window in [windows[0], windows[2]]:
            expected_ave += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"]) - (window["syn_subs"]/window["exp_syn_subs"]))*(window["nonsyn_subs"]+window["syn_subs"])/(window["treelen"]*exp_total_subs)

        actual_ave = site_dnds_info.get_ave_dn_minus_ds_weightby_subs(is_exclude_low_sub=True)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn-ds weighted by subs, exclude low sub sites=" + str(expected_ave) +
                                   " but got " + str(actual_ave))

        # Exclude low sub windows, weight dn-ds/treelen by reads
        exp_total_reads = sum([window["reads"] for window in [windows[0], windows[2]]])
        expected_ave = 0
        for window in [windows[0], windows[2]]:
            expected_ave += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"]) - (window["syn_subs"]/window["exp_syn_subs"]))*(window["reads"])/(window["treelen"]*exp_total_reads)

        actual_ave = site_dnds_info.get_ave_dn_minus_ds_weightby_reads(is_exclude_low_sub=True)
        self.assertAlmostEqual(expected_ave, actual_ave, places=10,
                               msg="Expected average dn-ds weighted by reads, exclude low sub sites=" + str(expected_ave) +
                                   " but got " + str(actual_ave))



    def test_site_dnds_info_(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows
        expected_windows = len(windows)
        actual_windows = site_dnds_info.get_window_coverage()
        self.assertEqual(expected_windows, actual_windows,
                               msg="Expected total windows covering site=" + str(expected_windows) +
                                   " but got " + str(actual_windows))

    def test_site_dnds_info_reads(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows
        expected_ave = sum([window["reads"] for window in windows])/len(windows)
        actual_ave = site_dnds_info.get_ave_read_coverage()
        self.assertEqual(expected_ave, actual_ave,
                               msg="Expected ave reads covering site=" + str(expected_ave) +
                                   " but got " + str(actual_ave))

    def test_site_dnds_info_nonsyn(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows
        expected_ave = sum([window["nonsyn"] for window in windows])/len(windows)
        actual_ave = site_dnds_info.get_ave_nonsyn_subs()
        self.assertEqual(expected_ave, actual_ave,
                               msg="Expected ave nonsynonymous substitutions at site=" + str(expected_ave) +
                                   " but got " + str(actual_ave))



    def test_site_dnds_info_syn(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows
        expected_ave = sum([window["syn"] for window in windows])/len(windows)
        actual_ave = site_dnds_info.get_ave_syn_subs()
        self.assertEqual(expected_ave, actual_ave,
                               msg="Expected ave synonymous substitutions at site=" + str(expected_ave) +
                                   " but got " + str(actual_ave))


    def test_site_dnds_info_sub(self):
        site_dnds_info = slice_miseq.SiteDnDsInfo()
        windows = self.windows
        expected_ave = sum([window["syn"]+window["nonsyn"] for window in windows])/len(windows)
        actual_ave = site_dnds_info.get_ave_subs()
        self.assertEqual(expected_ave, actual_ave,
                               msg="Expected ave substitutions at site=" + str(expected_ave) +
                                   " but got " + str(actual_ave))


    def test_tabulate_dnds(self):
        # Site-Window Coverage
        # Site  1   2   3   4
        # Win0  *   *
        # Win1      *   *   *
        # Win2      *   *
        # Win3      *   *
        tmpfilenames = []
        for window in self.windows:
            startnuc = (window["startcdn"]*3) - 2
            endnuc = window["endcdn"]*3
            tmpfile = tempfile.NamedTemporaryFile(delete=False, dir=OUT_DIR,
                                                  suffix=".{}_{}.dnds.tsv".format(startnuc, endnuc))
            writer = csv.DictWriter(tmpfile, fieldnames=[hyphy_handler.HYPHY_TSV_N_COL,
                                                         hyphy_handler.HYPHY_TSV_S_COL,
                                                         hyphy_handler.HYPHY_TSV_EXP_N_COL,
                                                         hyphy_handler.HYPHY_TSV_EXP_S_COL,
                                                         hyphy_handler.HYPHY_TSV_DN_COL,
                                                         hyphy_handler.HYPHY_TSV_DS_COL,
                                                         hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL],
                                    delimiter="\t")
            writer.writeheader()

            # Each window has the same values repeated for every site it covers
            for cdn_site in range(window["startcdn"], window["endcdn"]+1):
                outrow = dict()
                outrow[hyphy_handler.HYPHY_TSV_N_COL] = window["nonsyn_subs"]
                outrow[hyphy_handler.HYPHY_TSV_S_COL] = window["syn_subs"]
                outrow[hyphy_handler.HYPHY_TSV_EXP_N_COL] = window["exp_nonsyn_subs"]
                outrow[hyphy_handler.HYPHY_TSV_EXP_S_COL] = window["exp_syn_subs"]
                outrow[hyphy_handler.HYPHY_TSV_DN_COL] = window["nonsyn_subs"]/window["exp_nonsyn_subs"]
                outrow[hyphy_handler.HYPHY_TSV_DS_COL] = window["syn_subs"]/window["exp_syn_subs"]
                outrow[hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL] = ((window["nonsyn_subs"]/window["exp_nonsyn_subs"])-(window["syn_subs"]/window["exp_syn_subs"]))/window["treelen"]
                writer.writerow(outrow)

            tmpfile.flush()  # flush to os
            os.fsync(tmpfile.file.fileno())  # flush to disk
            tmpfile.close()

            # Create a fake window slice fasta.  We only care about the read depth is truthful.
            # We don't care if the fasta substitutions match the window substitutions
            tmpfastafile = tmpfile.name.replace(".dnds.tsv", ".fasta")
            with open(tmpfastafile, 'w') as fh_fasta:
                num_codons = window["endcdn"] - window["startcdn"] + 1
                for i in range(window["reads"]):
                    fh_fasta.write(">Read" + str(i) + "\n")
                    fh_fasta.write("ACT"*num_codons + "\n")

            tmpfilenames.extend([tmpfile.name, tmpfastafile])

        COMMENTS = "My comments"
        REF = "MyRef"
        REF_NUC_LEN = 12
        OUTPUT_CSV = OUT_DIR + os.sep + "test_tabulate_dnds.csv"
        slice_miseq.tabulate_dnds(dnds_tsv_dir=OUT_DIR, ref=REF,ref_nuc_len=REF_NUC_LEN,
                                  output_csv_filename=OUTPUT_CSV,
                                  comments=COMMENTS)

        # Check that the file exists
        self.assertTrue(os.path.exists(OUTPUT_CSV) and os.path.getsize(OUTPUT_CSV) > 0,
                        msg="Unable to write to dn/ds summary csv file")

        # verify contents
        with open(OUTPUT_CSV, 'rU') as fh_in_tabulate_csv:
            reader = csv.DictReader(filter(lambda row: row[0]!='#', fh_in_tabulate_csv))
            for site_0based, row in enumerate(reader):
                self.assertEqual(row["Site"], str(site_0based+1),
                                 "Expected Sites to be in order.  Expect site=" + str(site_0based+1) + " but got " + row["Site"])

                exp_total_subs_hisub = 0
                exp_total_subs_losub = 0
                exp_total_reads_hisub = 0
                exp_total_reads_losub = 0
                for winIndices in self.site_2_window[site_0based]:
                    if isinstance(winIndices, int ):
                        winIndices = [winIndices] # hack for when there is only 1 window covering a site

                    for windowIdx in winIndices:
                        window = self.windows[windowIdx]
                        if window["nonsyn_subs"] >= 1 and window["syn_subs"] >= 1:
                            exp_total_subs_hisub += window["nonsyn_subs"] + window["syn_subs"]
                            exp_total_reads_hisub += window["reads"]
                        else:
                            exp_total_subs_losub += window["nonsyn_subs"] + window["syn_subs"]
                            exp_total_reads_losub += window["reads"]

                expected_ave_dnds_weightByReads = 0
                expected_ave_dnds_weightByReads_nolowsub = 0
                expected_ave_dnds_weightBySubs = 0
                expected_ave_dnds_weightBySubs_nolowsub = 0

                expected_ave_dnminusds_weightByReads = 0
                expected_ave_dnminusds_weightByReads_nolowsub = 0
                expected_ave_dnminusds_weightBySubs = 0
                expected_ave_dnminusds_weightBySubs_nolowsub = 0
                for winIndices in self.site_2_window[site_0based]:
                    if isinstance(winIndices, int ):
                        winIndices = [winIndices] # hack for when there is only 1 window covering a site

                    for windowIdx in winIndices:
                        window = self.windows[windowIdx]

                        expected_ave_dnds_weightBySubs += window["nonsyn_subs"]*window["exp_syn_subs"]*(window["nonsyn_subs"]+window["syn_subs"])/(window["syn_subs"]*window["exp_nonsyn_subs"]*(exp_total_subs_hisub+exp_total_subs_losub))
                        expected_ave_dnds_weightByReads += window["nonsyn_subs"]*window["exp_syn_subs"]*(window["reads"])/(window["syn_subs"]*window["exp_nonsyn_subs"]*(exp_total_reads_hisub+exp_total_reads_losub))
                        expected_ave_dnminusds_weightByReads += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"])-(window["syn_subs"]/window["exp_syn_subs"]))*window["reads"]/(window["treelen"]*(exp_total_reads_hisub+exp_total_reads_losub))
                        expected_ave_dnminusds_weightBySubs += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"])-(window["syn_subs"]/window["exp_syn_subs"]))*(window["nonsyn_subs"]+window["syn_subs"])/(window["treelen"]*(exp_total_subs_hisub+exp_total_subs_losub))

                        if window["nonsyn_subs"] >= 1 and window["syn_subs"] >= 1:
                            expected_ave_dnds_weightBySubs_nolowsub += window["nonsyn_subs"]*window["exp_syn_subs"]*(window["nonsyn_subs"]+window["syn_subs"])/(window["syn_subs"]*window["exp_nonsyn_subs"]*exp_total_subs_hisub)
                            expected_ave_dnds_weightByReads_nolowsub += window["nonsyn_subs"]*window["exp_syn_subs"]*(window["reads"])/(window["syn_subs"]*window["exp_nonsyn_subs"]*exp_total_reads_hisub)
                            expected_ave_dnminusds_weightByReads_nolowsub += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"])-(window["syn_subs"]/window["exp_syn_subs"]))*window["reads"]/(window["treelen"]*exp_total_reads_hisub)
                            expected_ave_dnminusds_weightBySubs_nolowsub += ((window["nonsyn_subs"]/window["exp_nonsyn_subs"])-(window["syn_subs"]/window["exp_syn_subs"]))*(window["nonsyn_subs"]+window["syn_subs"])/(window["treelen"]*exp_total_subs_hisub)

                self.assertAlmostEqual(expected_ave_dnds_weightByReads, float(row["dndsWeightByReads"]), places=10,
                                       msg="Expected dn/ds weighted by reads=" + str(expected_ave_dnds_weightByReads) +
                                       " but got " + row["dndsWeightByReads"] + " at site " + str(site_0based+1))
                self.assertAlmostEqual(expected_ave_dnds_weightBySubs, float(row["dNdSWeightBySubs"]), places=10,
                                       msg="Expected dn/ds weighted by subs=" + str(expected_ave_dnds_weightBySubs) +
                                       " but got " + row["dnMinusDsWeightBySubs"] + " at site " + str(site_0based+1))
                self.assertAlmostEqual(expected_ave_dnminusds_weightByReads, float(row["dnMinusDsWeightByReads"]), places=10,
                                       msg="Expected dn-ds weighted by reads=" + str(expected_ave_dnminusds_weightByReads) +
                                       " but got " + row["dnMinusDsWeightByReads"] + " at site " + str(site_0based+1))
                self.assertAlmostEqual(expected_ave_dnminusds_weightBySubs, float(row["dnMinusDsWeightBySubs"]), places=10,
                                       msg="Expected dn-ds weighted by subs=" + str(expected_ave_dnminusds_weightBySubs) +
                                       " but got " + row["dnMinusDsWeightBySubs"] + " at site " + str(site_0based+1))

                if expected_ave_dnds_weightByReads_nolowsub:
                    self.assertAlmostEqual(expected_ave_dnds_weightByReads_nolowsub, float(row["dNdSWeightByReadsNoLowSub"]), places=10,
                                           msg="Expected dn/ds weighted by reads, excluding low sub sites=" + str(expected_ave_dnds_weightByReads_nolowsub) +
                                           " but got " + row["dNdSWeightByReadsNoLowSub"] + " at site " + str(site_0based+1))
                else:
                    self.assertEqual("", row["dNdSWeightByReadsNoLowSub"],
                                           msg="Expected dn/ds weighted by reads, excluding low sub sites=None"
                                           " but got " + row["dNdSWeightByReadsNoLowSub"] + " at site " + str(site_0based+1))
                if expected_ave_dnds_weightBySubs_nolowsub:
                    self.assertAlmostEqual(expected_ave_dnds_weightBySubs_nolowsub, float(row["dNdSWeightBySubsNoLowSub"]), places=10,
                                       msg="Expected dn/ds weighted by reads, excluding low sub sites=None" + str(expected_ave_dnds_weightBySubs_nolowsub) +
                                       " but got " + row["dNdSWeightBySubsNoLowSub"] + " at site " + str(site_0based+1))
                else:
                    self.assertEqual("", row["dNdSWeightBySubsNoLowSub"],
                                           msg="Expected dn/ds weighted by subs, excluding low sub sites=None"
                                           " but got " + row["dNdSWeightBySubsNoLowSub"] + " at site " + str(site_0based+1))

                if expected_ave_dnminusds_weightByReads_nolowsub:
                    self.assertAlmostEqual(expected_ave_dnminusds_weightByReads_nolowsub, float(row["dnMinusDsWeightByReadsNoLowSubs"]), places=10,
                                           msg="Expected dn-ds weighted by reads, excluding low sub sites=" + str(expected_ave_dnminusds_weightByReads_nolowsub) +
                                           " but got " + row["dnMinusDsWeightByReadsNoLowSubs"] + " at site " + str(site_0based+1))
                else:
                    self.assertEqual("", row["dnMinusDsWeightByReadsNoLowSubs"],
                                           msg="Expected dn-ds weighted by reads, excluding low sub sites=None"
                                           " but got " + row["dnMinusDsWeightByReadsNoLowSubs"] + " at site " + str(site_0based+1))
                if expected_ave_dnminusds_weightBySubs_nolowsub:
                    self.assertAlmostEqual(expected_ave_dnminusds_weightBySubs_nolowsub, float(row["dnMinusDsWeightBySubsNoLowSubs"]), places=10,
                                       msg="Expected dn-ds weighted by reads, excluding low sub sites=None" + str(expected_ave_dnminusds_weightBySubs_nolowsub) +
                                       " but got " + row["dnMinusDsWeightBySubsNoLowSubs"] + " at site " + str(site_0based+1))
                else:
                    self.assertEqual("", row["dnMinusDsWeightBySubsNoLowSubs"],
                                           msg="Expected dn-ds weighted by subs, excluding low sub sites=None"
                                           " but got " + row["dnMinusDsWeightBySubsNoLowSubs"] + " at site " + str(site_0based+1))




        for tmpfilename in tmpfilenames:
            os.remove(tmpfilename)

if __name__ == '__main__':
    unittest.main()
