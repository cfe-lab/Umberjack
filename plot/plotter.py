"""
Plots results along reference
"""
import subprocess
import os
import logging
import config.settings as settings


LOGGER = logging.getLogger(__name__)


RSCRIPT_EXE = "Rscript"

def plot_dnds(dnds_csv, dnds_img=None, Rscript_exe=RSCRIPT_EXE):
    """
    Uses R ggplot to plot the sitewise dN/dS averaged across windows to pdf.
    Tacks on .pdf to the csv
    :param str dnds_csv:  path to input per-site dN/dS CSV
    :param str dnds_img:  path to output pdf to store plot.
            If not supplied, then replaces the dnds_csv filename ".csv" with ".pdf" or tacks on ".pdf" to the dnds_csv filename.
    :param str Rscript_exe:  path to Rscript executable.  If empty, takes Rscript from PATH environment variable.
    :return:
    """
    LOGGER.debug("Plotting Ave dN/dS results for " + dnds_csv)
    if not dnds_img:
        if dnds_csv.endswith(".csv"):
            dnds_img = dnds_csv.replace(".csv", ".pdf")
        else:
            dnds_img = dnds_csv + ".pdf"
    rcmd = ["Rscript",
            os.path.dirname(os.path.realpath(__file__)) + os.sep + os.pardir + os.sep + "R" + os.sep + "plot_dNdS.R",
            dnds_csv,
            dnds_img]
    subprocess.check_call(rcmd, env=os.environ)
    LOGGER.debug("Done Plotting Ave dN/dS results for " + dnds_csv)
