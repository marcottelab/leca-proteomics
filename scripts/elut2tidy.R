# if(!require(argparse)){
#   install.packages("argparse")
#   suppressPackageStartupMessages(library(argparse))
# }

library(magrittr)
library(argparse)
library(dplyr)
library(tidyr)
library(readr)

#################################################
# THIS SCRIPT WAS DEVELOPED AND TESTED BY:
# Rachael M. Cox
# rachaelcox@utexas.edu
# Last updated: 6/7/2021
#################################################

# create parser object
parser <- ArgumentParser()
parser$add_argument("-i", "--elut_file", dest="elutfile", required=TRUE,
                    help="Input elutions file, .elut file output by msblender2elution.py")
parser$add_argument("-o", "--out_file", dest="outfile",
                    help="Outfile name for tidied (long format) elution data")
parser$add_argument("-s", "--sep", dest="sep", default=",", help="Separator for infile")
args <- parser$parse_args()

elut2tidy <- function(elutfile, outfile, sep){

  print("Input elution matrix:")
  elut_df <- read_delim(elutfile, delim = sep)
  print(head(elut_df, 10))
  
  tidy_df <- elut_df %>%
    rename(orthogroup = 1) %>%
    pivot_longer(!orthogroup, names_to = "fraction_id", 
                 values_to = "PSMs")
  
  print(paste("Tidy data output written to:", outfile))
  print(head(tidy_df, 25))
  
  write_csv(tidy_df, outfile)
  
}

elut2tidy(args$elutfile, args$outfile, args$sep)
