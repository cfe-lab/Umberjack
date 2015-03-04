# Imports configurations from /SlidingWindow/config/slidingwindow.config file.
from ConfigParser import SafeConfigParser
import os
import sys
import logging
import logging.config

# Logging Configs
###########################################
LOG_CONFIG_FILE = os.path.dirname(os.path.realpath(__file__)) + os.sep + "logging.conf"

# if the logging.conf file doesn't exist, then output all logging to stdout, set level to DEBUG
if not os.path.exists(LOG_CONFIG_FILE):
    ROOT_LOGGER = logging.getLogger("")
    ROOT_LOGGER.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s')
    console_handler.setFormatter(formatter)
    ROOT_LOGGER.addHandler(console_handler)

    ROOT_LOGGER.warn("Unable to find logging.conf file.  Using DEBUG level logging to stdout")
else:
    logging.config.fileConfig(LOG_CONFIG_FILE, disable_existing_loggers=False)

LOGGER = logging.getLogger(__name__)

# Execution Configs
###########################################
SECTION = "slidingwindow"
SLIDINGWINDOW_CONFIG_FILE = os.path.dirname(os.path.realpath(__file__)) + os.sep + "slidingwindow.config"

# Set defaults
if not os.path.exists(SLIDINGWINDOW_CONFIG_FILE):
    SAM_FILENAME = None
    REF = None
    OUT_DIR = None
    MAP_QUAL_CUTOFF = 20
    READ_QUAL_CUTOFF = 20
    MAX_PROP_N = 0.1
    WINDOW_SIZE = 300
    WINDOW_SLIDE = 30
    WINDOW_BREADTH_CUTOFF = 0.8
    WINDOW_DEPTH_CUTOFF = 50
    START_NUCPOS = 1
    END_NUCPOS = None
    THREADS_PER_WINDOW = 1
    CONCURRENT_WINDOWS = 1
    OUTPUT_CSV_FILENAME = "slidingwindow.dnds.csv"
    HYPHY_EXE = "HYPHYMP"
    HYPHY_BASEDIR = "/usr/local/lib/hyphy/TemplateBatchFiles/"
    FASTREE_EXE = "FastTreeMP"
    MODE =  "DNDS"
    MPI = False
    DEBUG = False
else:
    parser = SafeConfigParser(allow_no_value=True,
                              defaults={"sam_filename": None,
                                        "ref": None,
                                        "out_dir": None,
                                        "map_qual_cutoff":20,
                                        "read_qual_cutoff":20,
                                        "max_prop_n": 0.1,
                                        "window_size": 300,
                                        "window_slide": 30,
                                        "window_breadth_cutoff": 0.8,
                                        "window_depth_cutoff":50,
                                        "start_nucpos":1,
                                        "end_nucpos": None,
                                        "threads_per_window":1,
                                        "concurrent_windows":1,
                                        "mode": "DNDS",
                                        "hyphy_exe": "HYPHYMP",
                                        "hyphy_basedir": "/usr/local/lib/hyphy/TemplateBatchFiles/",
                                        "fastree_exe": "FastTreeMP",
                                        "output_csv_filename": "slidingwindow.out.csv",
                                        "mpi":  False,
                                        "debug": False,
                              })
    parser.readfp(fp=None, filename=SLIDINGWINDOW_CONFIG_FILE)

    SAM_FILENAME = parser.get(SECTION, "sam_filename")
    REF = parser.get(SECTION, "ref")
    OUT_DIR = parser.get(SECTION, "out_dir")
    MAP_QUAL_CUTOFF = parser.getint(SECTION, "map_qual_cutoff")
    READ_QUAL_CUTOFF = parser.getint(SECTION, "read_qual_cutoff")
    MAX_PROP_N = parser.getfloat(SECTION, "max_prop_n")
    WINDOW_SIZE = parser.getint(SECTION, "window_size")
    WINDOW_SLIDE = parser.getint(SECTION, "window_slide")
    WINDOW_BREADTH_CUTOFF = parser.getfloat(SECTION, "window_breadth_cutoff")
    WINDOW_DEPTH_CUTOFF = parser.getint(SECTION, "window_depth_cutoff")
    START_NUCPOS = parser.getint(SECTION, "start_nucpos")
    END_NUCPOS = parser.getint(SECTION, "end_nucpos")
    THREADS_PER_WINDOW = parser.getint(SECTION, "threads_per_window")
    CONCURRENT_WINDOWS = parser.getint(SECTION, "concurrent_windows")
    OUTPUT_CSV_FILENAME = parser.get(SECTION, "output_csv_filename")
    HYPHY_EXE = parser.get(SECTION, "hyphy_exe")
    HYPHY_BASEDIR = parser.get(SECTION, "hyphy_basedir")
    FASTREE_EXE = parser.get(SECTION, "fastree_exe")
    MODE = parser.get(SECTION, "mode")
    MPI = parser.getboolean(SECTION, "mpi")
    DEBUG = parser.getboolean(SECTION, "debug")

