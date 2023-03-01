#################################################
# THIS SCRIPT WAS DEVELOPED AND TESTED BY:
# Rachael M. Cox
# rachaelcox@utexas.edu
# Last updated: 01/18/2023
#################################################

library(argparse)
library(tidyverse)
parser <- ArgumentParser()
parser$add_argument("-i", "--cmplx_file", dest="cmplx_file", required=TRUE,
                    help="File containing a complex (one row per protein or orthogroup), identifier column = 'ID'")
parser$add_argument("-o", "--outfile", dest="outfile",
                    help="Outfile path/name for dot plot figures (not required)")
parser$add_argument("-s", "--sep", dest="sep", default=",", help="Separator for infile (default = ',')")
args <- parser$parse_args()

# -----------------------------------------------------
################# global data set up ##################
# -----------------------------------------------------

# coverage data
cov <- read_csv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/coverage/og_coverage_by_species_phylo.150p.filtdollo.csv')

# species names & clades
species <- read_csv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/speciesinfo_clades.csv') %>%
  filter(tax_group == "eukaryota") %>% 
  select(code, species_name, clade) %>%
  mutate(code = tolower(code))

# ordered codes based on phylogeny
species_order <- read_delim('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/euk_codes_ordered_phylo.txt',
                   delim='\n', col_names=FALSE) %>% pull(X1)

# clade order
clade_order <- c('Amorphea','Excavate','TSAR','Archaeplastida')

# og annotations
annots <- read_tsv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/leca_euNOGs_annotated.fmt.tsv') %>%
  select(ID, matches('*genes_fmt*'),
         human_protein_names,
         human_function_cc) %>%
  mutate(gene_names = coalesce(human_genes_fmt,
                               arath_genes_fmt,
                               ID)) %>%
  select(-matches('*genes_fmt*')) %>%
  select(ID, gene_names, everything())

multimapped <- annots %>%
  group_by(gene_names) %>%
  tally() %>%
  filter(n > 1)

annots$gene_names = gsub(x = annots$gene_names, pattern = ";",
                             replacement = ",")
annots$gene_names = gsub(x = annots$gene_names, pattern = "NA ",
                             replacement = "")

annots_fmt <- annots %>%
  mutate(gene_names = str_replace(gene_names, "^([^,]*,[^,]*),.*", "\\1 ...")) %>% 
  mutate(id_suffix = paste0("(",ID,")")) %>%
  mutate(gene_names = ifelse(gene_names %in% multimapped$gene_names, 
                             paste(gene_names,id_suffix,sep="\n"), 
                             gene_names)) %>%
  select(-id_suffix)

# -----------------------------------------------------
#################### functions ####################
# -----------------------------------------------------

# get coverage data for complex
get_phylo_data <- function(cmplx_file, sep = ','){
  
  # read in desired complx
  cmplx <- read_delim(cmplx_file, delim = sep)
  
  cmplx_name <- cmplx %>%
    pull(granulated_cmplx_name) %>%
    unique()
  
  # get cluster order for prots
  cmplx_ids_ordered <- cmplx %>%
    select(ID, characterization_status) %>%
    unique() %>%
    mutate(fct_lvl = row_number())
  
  message("Pulling phylogenetic information for: ", cmplx_name)
  cov_cmplx <- cov %>%
    filter(ID %in% cmplx$ID) %>%
    left_join(cmplx_ids_ordered)
  
  print(cov_cmplx)
  return(cov_cmplx)
  
}

# format data for plotting
fmt_phylo_data <- function(cov_cmplx){
  
  message("\nFormatting data ...")
  # tidy data & join annots
  cov_tidy <- cov_cmplx %>%
    pivot_longer(-c(ID, characterization_status, fct_lvl),
                 names_to = "species", values_to = "presence") %>%
    left_join(species, by=c("species" = "code")) %>%
    left_join(annots_fmt, by=c("ID")) %>%
    mutate(gene_names_nov = ifelse(str_detect(characterization_status, 'Novel'),
                                   paste0('*', gene_names),
                                   gene_names)) %>%
    mutate(gene_names = fct_reorder(gene_names, fct_lvl)) %>% 
    mutate(gene_names_nov = fct_reorder(gene_names_nov, fct_lvl)) # maintain prot cluster order
  
  # lock in other var orders
  cov_tidy$clade <- factor(cov_tidy$clade, levels = clade_order)
  cov_tidy$species <- factor(cov_tidy$species, levels = species_order)
  
  return(cov_tidy)
  
}

# generate dot plot
plot_dots <- function(cov_tidy){
  
  message("\nGenerating plot ...")
  
  pal_clades <- c("#E64B35", "#4DBBD5", "#3C5488", "#00A087")
  
  # get cmplx size
  num_units <- length(pull(cov_tidy, ID) %>% unique)
  
  # adjust dots based on cmplx size
  if(num_units > 12){
    dot_size = 1
  } else {
    dot_size = 2.5
  }
  
  p <- ggplot(cov_tidy, aes(x = species, y = gene_names_nov, color = clade)) +
    geom_point(shape = 1, color = "black") +
    geom_point(aes(color = clade, 
                   size = ifelse(presence==0, NA, dot_size))) +
    theme_bw() +
    scale_color_manual(values = pal_clades) +
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.margin=margin(t=-10)) +
    scale_size(guide = 'none') +
    labs(color = "")
  
  return(p)
  
}

# -----------------------------------------------------
#################### wrapper ####################
# -----------------------------------------------------

plot_phylo_dots <- function(cmplx_file, sep = ',',
                            outfile = NULL){
  
  phylo_data <- get_phylo_data(cmplx_file)
  phylo_fmt <- fmt_phylo_data(phylo_data)
  plot <- plot_dots(phylo_fmt)
  
  if(!is.null(outfile)){
    outfile_name <- outfile
  } else {
    outfile_name <- tools::file_path_sans_ext(cmplx_file)
    outfile_name <- paste0(outfile_name, "_phyloplot")
  }
  
  message("Saving plots to:\n", 
          paste0(outfile_name, ".pdf\n"), 
          paste0(outfile_name, ".png"))
  
  # get cmplx size
  num_units <- length(pull(phylo_data, ID) %>% unique)
  print(paste0("# subunits = ", num_units))
  
  # save plot dimensions based on cmplx size
  if(num_units >= 12){
    
    print(paste0("# subunits = ", num_units))
    message("Adjusting plot dimensions ...")
    hval = 0.33*num_units
    hval_max = 24
    
    plot %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                    width = 6, height = ifelse(hval < hval_max, hval, hval_max),
                    units = "in")
    plot %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                    width = 6, height = ifelse(hval < hval_max, hval, hval_max),
                    units = "in")
    
  } else {
    
    plot %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                    width = 6, height = 3, units = "in")
    plot %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                    width = 6, height = 3, units = "in")
    
    message("Done!")
    
  }
  
  plot %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                  width = 6, height = 3, units = "in", dpi = 300)
  plot %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                  width = 6, height = 3, units = "in", dpi = 300)
  
  message("Done!")
  
}

plot_phylo_dots(cmplx_file = args$cmplx_file, outfile = args$outfile, sep = args$sep)