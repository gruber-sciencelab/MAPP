###############################################################################
#
#   Create .png files with sequence logos for provided motifs.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 21-03-2020
#   LICENSE: Apache_2.0
#   USAGE:
#   Rscript create-seqlogos.r --PWM_dir {dir} --output_dir {dir}
#
###############################################################################

# by default: suppress warnings
options(warn = -1)

# load libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("seqLogo"))

# list the command-line arguments
option_list <- list(
  make_option(c("--PWM_dir"),
    action = "store",
    dest = "PWM_dir",
    type = "character",
    default = ".",
    help = "Path to the directory with PWM files."
  ),

  make_option(c("--output_dir"),
    action = "store",
    dest = "output_dir",
    type = "character",
    default = "seqlogos",
    help = "Path to the directory for the sequence logos."
  ),

  make_option(c("--help"),
    action = "store_true",
    dest = "help",
    type = "logical",
    default = FALSE,
    help = "Show this information and exit."
  ),

  make_option(c("--verbose"),
    action = "store_true",
    dest = "verbose",
    type = "logical",
    default = FALSE,
    help = "Run in verbose mode."
  )
)

# parse command-line arguments
opt_parser <- OptionParser(
  usage = "Usage: %prog --PWM_dir [PATH] --output_dir [PATH]",
  option_list = option_list,
  add_help_option = FALSE,
  description = ""
)
opt <- parse_args(opt_parser)

# if verbose flag was set: print warnings
if (opt$verbose) {
  options(warn = 0)
}

###############################################################################
# MAIN
###############################################################################

PWM_paths <- Sys.glob(file.path(opt$PWM_dir, "*"))

for (path in PWM_paths) {
  fname <- basename(path)
  outfile <- paste(opt$output_dir, paste(fname, "png", sep = "."), sep = "/")

  # nice hack to read MetEvo PWMs format as R data frame
  pwm_df <- read.table(textConnection(rev(rev(readLines(path))[-1])), skip = 2)
  pwm_df <- t(pwm_df / 100.0)
  pwm_pwm <- makePWM(pwm_df, alphabet = "RNA")

  png(filename = outfile, height = 500, width = 800, pointsize = 10)
  seqLogo(pwm_pwm)
  dev.off()
}
