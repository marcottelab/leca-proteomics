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
parser <- argparse::ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-d", "--file_dir", dest="file_dir", required=TRUE,
                    help="Path to directory containing Complex Portal files that are tab-separated")
parser$add_argument("-p", "--file_pattern", dest="file_pattern", required=TRUE,
                    help="Pattern of files to format (e.g., '*tsv')")
parser$add_argument("-a", "--write_annotations", dest="write_annot", default=FALSE,
                    help="Choice for outputting annotation files in addition to formatted complexes; default=FALSE")

args <- parser$parse_args()

# function for formatting complex IDs and writing out gold std protein interaction file
fmt_cmplx <- function(cmplx_prtl_file, outfile_suffix, write_annot = FALSE){
  
  if(missing(outfile_suffix)){
    
    outfile_suffix = ".gold.cmplx.txt"
    
  }
  
  filename <- tools::file_path_sans_ext(cmplx_prtl_file)
  outfilename <- paste0(filename, outfile_suffix)
  
  cmplx_portal <- read_tsv(cmplx_prtl_file) %>%
    janitor::clean_names()
  
  cmplx_fmt <- cmplx_portal %>%
    # select(number_complex_ac, recommended_name,
    #        identifiers_and_stoichiometry_of_molecules_in_complex,
    #        expanded_participant_list, complex_assembly) %>%
    rename(identifiers = identifiers_and_stoichiometry_of_molecules_in_complex) %>%
    mutate(ids_fmt = str_replace_all(expanded_participant_list, "\\|", " ")) %>%
    mutate(ids_fmt = str_replace_all(ids_fmt, "\\(\\d\\)", "")) %>%
    mutate(ids_fmt = str_replace_all(ids_fmt, "-PRO.*?\\b", "")) %>% 
    mutate(ids_fmt = str_replace_all(ids_fmt, "-\\d", "")) %>%
    select(ids_fmt, everything())
  
  if(write_annot == TRUE){
    
    outfilename_annot <- paste0(filename, ".annotations", outfile_suffix)
    write_csv(cmplx_fmt, outfilename_annot, col_names = TRUE)
  
    }
    
  cmplx_final <- select(cmplx_fmt, ids_fmt)
  write_delim(cmplx_final, outfilename, delim = "\n", col_names = FALSE)
  
}

# locate all data
print("Locating data...")
file_path = args$file_dir
file_pattern = args$file_pattern
write_annot = args$write_annot

sprintf("Looking for Complex Portal files in %s/%s ...", file_path, file_pattern)

cmplx_files <- dir(file_path, pattern = file_pattern, full.names = TRUE)
sprintf("Processing the following %f file(s) ...", length(cmplx_files))
print(cmplx_files)

# format Complex Portal files
cmplxes_fmt <- mapply(fmt_cmplx,
                      cmplx_prtl_file = cmplx_files,
                      write_annot = write_annot)