suppressPackageStartupMessages(library(tidyverse))

if(!require(argparse)){
  install.packages("argparse")
  suppressPackageStartupMessages(library(argparse))
}

parser <- ArgumentParser()

parser$add_argument("-f", "--file", dest="file", required=TRUE,
                    help="Input file with UniProt accessions that need to be formatted; e.g., 'sp|A0A024RBG1|NUD4B_HUMAN' -> 'A0A024RBG1'. Column should be named 'ProteinID'.")

args <- parser$parse_args()

read_tsv(args$file) %>%
  mutate(ProteinID = str_extract(ProteinID,'(?<=\\|)(.*)(?=\\|)')) %>%
  write_tsv(paste0(args$file, ".fmt")) %>%
  print(args$file)