#!/usr/bin/env Rscript

if(!require(argparse)){
  install.packages("argparse")
  suppressPackageStartupMessages(library(argparse))
}

if(!require(tidyverse)){
  install.packages("tidyverse")
  suppressPackageStartupMessages(library(tidyverse))
}

# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-f", "--featmat", dest="featmat", required=TRUE,
                    help="Comma-separated feature matrix where the first two columns correspond to ID1, ID2")
# object for parsing input variables
args <- parser$parse_args()

featmat = args$featmat
filename <- tools::file_path_sans_ext(featmat)
outfile_suffix = ".summarized.featmat"
outfilename <- paste0(filename, outfile_suffix)

group_ids <- function(featfile, outfile){
  print("Reading in feature matrix ... ")
  featmat_dt <- data.table::fread(featfile, sep = ",")  # read in feature matrix
  featmat_dt[is.na(featmat_dt)] <- 0  # convert missing values to 0
  print("Grouping by ID1, ID2 and calculating the mean of each feature ...")
  summarized <- featmat_dt %>%
    group_by(acc1, acc2) %>%
    summarise(across(everything(), mean))
  
  print(paste0("Writing results to: ", outfile))
  write_csv(summarized, outfile)
}

group_ids(featmat, outfilename)