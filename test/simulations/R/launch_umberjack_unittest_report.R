# We can't pass arguments directly to an R script that we want to spin to a knitr HTML doc.
# But the spun R script will persist existing R environment variables.
library(knitr)
args <- commandArgs(TRUE)
if (length(args) < 2) stop(paste0("Bad args, usage:", 
                                  "Rscript  'launch_umberjack_unittest_report.R' [filepath to report config file] [filepath of output HTML report]"))

CONFIG_FILENAME <- args[1]
REPORT_OUT_DIR <- args[2]

# Set the directory for all the output figures
opts_knit$set(base.dir = REPORT_OUT_DIR)
# Must create the figure subdirectory if we use a custom output directory
dir.create(paste0(REPORT_OUT_DIR, .Platform$file.sep, "figure"), recursive = TRUE, showWarnings = FALSE)

spin("umberjack_unit_test.R", knit=FALSE)
knit2html("umberjack_unit_test.Rmd", output=paste0(REPORT_OUT_DIR, .Platform$file.sep, "umberjack_unit_test.html"), stylesheet="markdown_bigwidth.css")
