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
parser$add_argument("-d", "--elut_dir", dest="elut_dir", required=TRUE,
                    help="Path to directory containing elut files that are tab-separated")
parser$add_argument("-p", "--elut_pattern", dest="elut_pattern",
                    help="Pattern for suffix of .elut files to normalize (e.g. *unique.elut)")

args <- parser$parse_args()

# function for reading, parsing and normalizing elut files
parse_elut <- function(elutfile, norm = FALSE, write = FALSE, outfile_suffix = ".parse.elut"){
  
  elut_df <- read_delim(elutfile, delim = "\t")
  print(elut_df)
  print(elut_df[1])
  
  if(norm == TRUE){
    
    elut_matrix <- elut_df %>% 
      select(-1)  # remove identifier column
    
    elut_norm <- 
      t(apply(elut_matrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
    
    elut_out <- data.frame(elut_norm, orthogroup = elut_df[1]) # add back on identifiers
    
    elut_df <- elut_out %>%
      select(orthogroup, everything())
  }
  
  if(write == TRUE){
    
    filename <- tools::file_path_sans_ext(elutfile)
    outfilename <- paste0(filename, outfile_suffix)
    write_csv(elut_df, outfilename)
    
  }
  
  return(elut_df)
}

# locate all data
elut_path = args$elut_dir
elut_pattern = args$elut_pattern

elut_files <- dir(elut_path, pattern = elut_pattern,
                  full.names = TRUE)
print("Processing the following files:")
print(elut_files)

# normalize each elution profile by experiment
results_norm <- mapply(parse_elut, 
                       elutfile = elut_files,
                       norm = TRUE,
                       write = TRUE,
                       outfile_suffix = ".norm.elut")
