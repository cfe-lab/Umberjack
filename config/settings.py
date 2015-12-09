# Imports configurations from /Umberjack/config/umberjack.config file.
from ConfigParser import SafeConfigParser
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import chain
import os
import logging
import logging.config
import sys

VERSION = "1.0.0"

# Logging Configs
###########################################
DEFAULT_LOG_CONFIG_FILE = os.path.dirname(os.path.realpath(__file__)) + os.sep + "logging.conf"



DEFAULT_LOG_CONFIG_DICT =  {
    "version": 1,
    "disable_existing_loggers": False,
    "filters": {
        "hostnamefilter": {
            "()": "config.filters.HostnameFilter"
        }
    },
    "formatters":{
        "full":{
            "format":"%(asctime)s - [%(levelname)s] [%(name)s] [%(hostname)s] [%(process)d] %(message)s",
            # "format":"%(asctime)s - [%(levelname)s] [%(name)s] [%(process)d] %(message)s",
            "datefmt": ""  # leave datefmt emtpy to use default ISO8601 format %Y-%m-%d %H:%M:%S,%s
        }
    },
    "handlers":{
        "rotatingfile": {
            "formatter": "full",
            "filters": ["hostnamefilter"],
            "class": "logging.handlers.RotatingFileHandler",
            "filename": "umberjack.log",
            "maxBytes": 1000000000,  # 1 GB
            "backupCount": 20,
            "mode": "a",
            "encoding": "utf8"
        },
        "console": {
            "formatter": "full",
            "filters": ["hostnamefilter"],
            "class": "logging.StreamHandler",
            "stream": sys.stdout,
        }
    },
    # In dictConfig, the root logger must be specified parallel to handlers, formatters, filters.
    # If the root logger config is inside the "loggers" tag, it will not apply to the root logger.
    "root":{
        "handlers":["rotatingfile", "console"],
        "level":"DEBUG"
    }
}


# Execution Configs
###########################################
UMBERJACK_CONFIG_SECTION = "umberjack"
DEFAULT_UMBERJACK_CONFIG_FILE = os.path.dirname(os.path.realpath(__file__)) + os.sep + "umberjack.config"


DEFAULT_SAM_FILENAME = None
DEFAULT_REF = None
DEFAULT_OUT_DIR = os.getcwd()
DEFAULT_MAP_QUAL_CUTOFF = 20
DEFAULT_READ_QUAL_CUTOFF = 20
DEFAULT_MAX_PROP_N = 0.1
DEFAULT_WINDOW_SIZE = 300
DEFAULT_WINDOW_SLIDE = 30
DEFAULT_WINDOW_BREADTH_CUTOFF = 0.8
DEFAULT_WINDOW_DEPTH_CUTOFF = 50
DEFAULT_START_NUCPOS = 1
DEFAULT_END_NUCPOS = 0
DEFAULT_THREADS_PER_WINDOW = 1
DEFAULT_CONCURRENT_WINDOWS = 1
DEFAULT_OUTPUT_CSV_FILENAME = "umberjack.out.csv"
DEFAULT_HYPHY_EXE = "HYPHYMP"
DEFAULT_HYPHY_BASEDIR = os.path.abspath(os.path.realpath(__file__) + os.sep + os.pardir + os.sep + os.pardir + os.sep + "hyphy" + os.sep + "batchfile")
DEFAULT_HYPHY_LIBDIR = os.path.abspath(os.path.realpath(__file__) + os.sep + os.pardir + os.sep + os.pardir + os.sep + "test" + os.sep + "simulations" + os.sep + "bin" + os.sep + "hyphy" + os.sep + "hyphy_2.2.3" + os.sep + "res")
DEFAULT_FASTREE_EXE = "FastTree"
DEFAULT_FASTTREEMP_EXE = "FastTreeMP"
DEFAULT_MODE =  "DNDS"
DEFAULT_MPI = False
DEFAULT_DEBUG = False
DEFAULT_MASK_STOP_CODON = True
DEFAULT_REMOVE_DUPLICATES = True
DEFAULT_INSERT = False


def setup_logging(config_file=DEFAULT_LOG_CONFIG_FILE, config_dict=DEFAULT_LOG_CONFIG_DICT):
    """
    Set up logging globally using config file.
    The priority for configurations:
    1) config_file
    2) If config file is not specified, then DEFAULT_LOG_CONFIG_FILE is used as config file
    3) If the config file is explicity set to None/empty string, then config_dict is used
    4) If config_dict is not specified, then DEFAULT_LOG_CONFIG_DICT is used as config_dict
    5) If config_file and config_dict are both set to None or empty, then no configurations are set for logging

    :param str config_file:  file path to logging configuration file.
    :param dict config_dict:  dict configuration for logging
    """
    if config_file:
        logging.config.fileConfig(config_file, disable_existing_loggers=False)

    elif config_dict:
        logging.config.dictConfig(config_dict)



def setup_umberjack_config(config_file=DEFAULT_UMBERJACK_CONFIG_FILE, argname_prefix=""):
    """
    Parse the umberjack config file to obtain equivalent commandline arguments.
    I.e.  converts ConfigParser options into a list of commandline arguments that ArgumentParser can understand.
    :param str config_file: file path to umberjack config file that contains program arguments.
    :param str argname_prefix:  prefix to prepend to the program argument keys
    :return list: list of equivalent commandline arguments
    """
    default_config_dict = {
        "sam_filename": DEFAULT_SAM_FILENAME,
        "ref": DEFAULT_REF,
        "out_dir": DEFAULT_OUT_DIR,
        "map_qual_cutoff":DEFAULT_MAP_QUAL_CUTOFF,
        "read_qual_cutoff":DEFAULT_READ_QUAL_CUTOFF,
        "max_prop_n": DEFAULT_MAX_PROP_N,
        "window_size": DEFAULT_WINDOW_SIZE,
        "window_slide": DEFAULT_WINDOW_SLIDE,
        "window_breadth_cutoff": DEFAULT_WINDOW_BREADTH_CUTOFF,
        "window_depth_cutoff": DEFAULT_WINDOW_DEPTH_CUTOFF,
        "start_nucpos": DEFAULT_START_NUCPOS,
        "end_nucpos": DEFAULT_END_NUCPOS,
        "insert": DEFAULT_INSERT,
        "mask_stop_codon": DEFAULT_MASK_STOP_CODON,
        "remove_duplicates": DEFAULT_REMOVE_DUPLICATES,
        "threads_per_window": DEFAULT_THREADS_PER_WINDOW,
        "concurrent_windows": DEFAULT_CONCURRENT_WINDOWS,
        "mode": DEFAULT_MODE,
        "hyphy_exe": DEFAULT_HYPHY_EXE,
        "hyphy_basedir": DEFAULT_HYPHY_BASEDIR,
        "hyphy_libdir": DEFAULT_HYPHY_LIBDIR,
        "fastree_exe": DEFAULT_FASTREE_EXE,
        "output_csv_filename": DEFAULT_OUTPUT_CSV_FILENAME,
        "mpi":  DEFAULT_MPI,
        "debug": DEFAULT_DEBUG
    }
    arg_list = []
    if os.path.exists(config_file):
        config_parser = SafeConfigParser(allow_no_value=True, defaults=default_config_dict)
        config_parser.read(config_file)
        for option in config_parser.options(UMBERJACK_CONFIG_SECTION):
            print str(option)
            arg_list.extend([argname_prefix + option, config_parser.get(UMBERJACK_CONFIG_SECTION, option)])
    else:
        for key, val in default_config_dict:
            arg_list.extend([argname_prefix + key, val])
    return arg_list

