"""
Handles the metadata for a simulated dataset made via sim_pipeline.py
"""
import re
import glob
import os
from ConfigParser import RawConfigParser

class SimData:

    SECTION = "sim"

    def __init__(self, config_file):
        self.config_file = config_file
        self.sim_data_dir = os.path.abspath(os.path.dirname(config_file))

        config = RawConfigParser()
        config.read(config_file)

        self.name = config.get(SimData.SECTION, "FILENAME_PREFIX")


    @staticmethod
    def get_recombo_breaks_from_file(filename):
        """
        Returns the recombinant section start and end in 1-based nucleotide coordinates.
        :param str filename: filepath to any file for a recombinant section in the sim_pipeline.py dataset generation.
        :return int, int:  section start, section end or (None, None) if no match found
        """
        # *.break.<section nuc 1based start>_<section nuc 1based end>.*
        section_start = None
        section_end = None
        section_match = re.search(pattern=r"\.break\.(\d+)_(\d+)\.", string=filename)
        if section_match:
            section_start = int(section_match.group(1))
            section_end = int(section_match.group(2))
        return section_start, section_end


    def get_recombo_trees(self):
        """
        Gets the newick treefiles whose branch lengths have been inferred from the full population sequence from FastTree
        :return list: list of str filepaths to breakpoint trees, sorted by breakpoint order
        """
        if not os.path.exists(self.sim_data_dir + os.sep + "subs"):
            raise ValueError("Must execute sim_pipeline.py.  The dataset is not yet done.")

        # There might be recombination in the full population, in which there will be multiple trees.
        full_popn_treefiles = glob.glob(self.sim_data_dir + os.sep + "subs" + os.sep + "*.break.*_*.fasttree.nwk")

        full_popn_treefiles = sorted(full_popn_treefiles, key=lambda treefile: SimData.get_recombo_breaks_from_file(treefile)[0])

        return full_popn_treefiles


    def get_recombo_breaks(self):
        """
        Returns the recombinant sections ranges
        :return list:  [(int section nuc start 1based, int section nuc end 1based)]
        """
        if not os.path.exists(self.sim_data_dir + os.sep + "subs"):
            raise ValueError("Must execute sim_pipeline.py.  The dataset is not yet done.")

        # There might be recombination in the full population, in which there will be multiple trees.
        section_ranges = []
        for full_popn_treefile in  glob.glob(self.sim_data_dir + os.sep + "subs" + os.sep + "*.break.*_*.fasttree.nwk"):
            section_ranges.append(SimData.get_recombo_breaks_from_file(full_popn_treefile))

        return section_ranges



    def get_fasta(self):
        """
        Gets the multiple sequence aligned full population fasta for the entire genome
        :return str:  filepath
        """
        if not os.path.exists(self.sim_data_dir + os.sep + "fullpopn"):
            raise ValueError("Must execute sim_pipeline.py.  The dataset is not yet done.")

        return self.sim_data_dir + os.sep + "fullpopn" + os.sep + self.name + "_TRUE.fasta"