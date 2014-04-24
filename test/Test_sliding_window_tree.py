import unittest
import sliding_window_tree
import sys, os
import array
import csv
import subprocess

# For now, use a *.remap.sam file (paired end reads aligned to a consensus sequence with indels removed).
# SAM_FILENAME = "./data/TestSample-RT_S17.HIV1B-vif.remap.sam"
# MAPQ_CUTOFF = 0  # alignment quality cutoff
# MAX_PROP_N = 0  # maximum proportion of N bases in MSA-aligned sequence
# READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
# REFERENCE_FASTA = "./data/TestSample-RT_S17.HIV1B-vif.10.conseq"

SAM_FILENAME = "./simulations/data/sample_genomes.grinder-reads.rename.sam"
MAPQ_CUTOFF = 0  # alignment quality cutoff
MAX_PROP_N = 0.2  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]
REFERENCE_FASTA = "./simulations/data/sample_genomes.consensus.fas"
OUT_DIR = "./simulations/data/out"

GAMMA_DNDS_LOOKUP_FILENAME = "./simulations/data/sample_genomes.rates"
GAMMA_DNDS_LOOKUP_COL_INTERVAL = "Interval"
GAMMA_DNDS_LOOKUP_COL_SITE = "Site"
GAMMA_DNDS_LOOKUP_COL_SCALE_FACTOR = "Scaling_factor"
GAMMA_DNDS_LOOKUP_COL_RATE_CLASS = "Rate_class"
GAMMA_DNDS_LOOKUP_COL_DNDS = "Omega"

INDELIBLE_RATE_OUTPUT_COLS = ["Site", "Class", "Partition", "Inserted?"]
INDELIBLE_RATE_OUTPUT_FILENAME_PATTERN = './simulations/data/scaling_{}.0_RATES.txt'
TOTAL_REF_CODON_SITES = 10000

SCALING_FACTORS = [1, 2, 5, 10, 20, 50, 100]


MIN_WINDOW_BREADTH_COV_FRACTION = 0.5
MIN_WINDOW_DEPTH_COV = 2

BWA_EXE = "/home/thuy/programs/bwa/bwa-0.7.8/bwa"  # TODO:  do not hardcode
BOWTIE_EXE = ""
BOWTIE_BUILD_EXE=""



class TestSlidingWindowTree(unittest.TestCase):

    def __get_expected_dnds(self, dnds_lookup_filename):
        """
        Look dn/ds from sample_genomes.rates file.
        :rtype array of float: array where each element is a site dn/ds value.
        :param indelible_rates_filename str: full file path to indelible rates output file.
        :param scaling.factor int: the amount the the branch lengths are scaled
        """

        with open(dnds_lookup_filename, 'r') as lookup_fh:
            site0based_to_dnds = array.array('f')
            for row in csv.DictReader(lookup_fh, delimiter=','):
                site = row[GAMMA_DNDS_LOOKUP_COL_SITE]
                dnds = row[GAMMA_DNDS_LOOKUP_COL_DNDS]
                site0based_to_dnds.append(float(dnds))

            return site0based_to_dnds


    # def __align(self, ref_fasta_filename, logfilename, threads, query_fastq_filename):
    #     query_filename_prefix = os.path.splitext(query_fastq_filename)[0]
    #     sai_filename = query_filename_prefix + ".sai"   # bwa aln output file
    #     with open(logfilename, 'w') as logfile_fh, open(sai_filename, 'w') as sai_fh:
    #         subprocess.check_call([BWA_EXE, 'index', ref_fasta_filename], stderr=logfile_fh, stdout=logfile_fh)
    #         subprocess.check_call([BWA_EXE, 'aln', '-t', threads, ref_fasta_filename, query_fastq_filename],
    #                               stderr=logfile_fh, stdout=sai_fh)
    #
    #         # TODO:  call this for pair 1 and pair 2
    #         subprocess.check_call([BWA_EXE, 'aln', '-t', threads, ref_fasta_filename, query_fastq_filename],
    #                               stderr=logfile_fh, stdout=sai_fh)
    #
    #         subprocess.check_call([BWA_EXE, 'sampe', '-P', ref_fasta_filename, ])
    #         sampe 	bwa sampe [-a maxInsSize] [-o maxOcc] [-n maxHitPaired] [-N maxHitDis] [-P] <in.db.fasta> <in1.sai> <in2.sai> <in1.fq> <in2.fq> > <out.sam>

    def __align(self, ref_fasta_filename, logfilename, threads, query_fq_1_filename, query_fq_2_filename):
        MEAN_FRAG_LEN = 500
        MAX_FRAG_LEN = MEAN_FRAG_LEN * 0,25
        sam_filename_prefix = os.path.splitext(logfilename)[0]

        subprocess.check_call(BOWTIE_BUILD_EXE, ref_fasta_filename)
        subprocess.check_call(BOWTIE_EXE, '--local', '-p', threads,
                              '-X', MAX_FRAG_LEN,
                              '-x', ref_fasta_filename,
                              '-1', query_fq_1_filename,
                              '-2', query_fq_2_filename,
                              '-S', sam_filename_prefix)


    def test_process_windows(self):
        
        actual_dnds_filename = './simulations/data/actual_dnds.tsv'
        with open(actual_dnds_filename, 'w') as dnds_fh:
            dnds_fh.write("ref\tsite\tdnds")
            expected_site_2_dnds = self.__get_expected_dnds(GAMMA_DNDS_LOOKUP_FILENAME)

            ref2SeqDnDsInfo = sliding_window_tree.process_windows(sam_filename=SAM_FILENAME,
                                                                ref_fasta_filename=REFERENCE_FASTA,
                                                                out_dir=OUT_DIR,
                                                                mapping_cutoff=MAPQ_CUTOFF,
                                                                read_qual_cutoff=READ_QUAL_CUTOFF,
                                                                max_prop_N=MAX_PROP_N,
                                                                window_breadth_thresh=MIN_WINDOW_BREADTH_COV_FRACTION,
                                                                window_depth_thresh=MIN_WINDOW_DEPTH_COV)

            for ref in ref2SeqDnDsInfo:
                for site in range(ref2SeqDnDsInfo[ref].get_seq_len()):
                    site_dnds = ref2SeqDnDsInfo[ref].get_site_ave_dnds(site)
                    dnds_fh.write(ref + "\t" + site + "\t" + site_dnds + "\n")



            
            

if __name__ == '__main__':
    unittest.main()
