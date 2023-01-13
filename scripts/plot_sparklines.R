#################################################
# THIS SCRIPT WAS DEVELOPED AND TESTED BY:
# Rachael M. Cox
# rachaelcox@utexas.edu
# Last updated: 01/11/2023
#################################################
library(argparse)
library(tidyverse)
parser <- ArgumentParser()
parser$add_argument("-i", "--cmplx_file", dest="cmplx_file", required=TRUE,
                    help="File containing a complex (one row per protein or orthogroup), identifier column = 'ID'")
parser$add_argument("-o", "--outfile", dest="outfile",
                    help="Outfile path/name for sparkline figures (not required)")
parser$add_argument("-s", "--sep", dest="sep", default=",", help="Separator for infile (default = ',')")
parser$add_argument("-e", "--best_exp", dest="best_exp", default=TRUE, help="Sub-sample each species for the best fractionation for the given complex (default = TRUE)")
parser$add_argument("-c", "--sample_clades", dest="sample_clades", default=FALSE, help="Sub-sample each eukaryotic clade for the given complex (default = FALSE)")
args <- parser$parse_args()

# -----------------------------------------------------
################# global data set up ##################
# -----------------------------------------------------
# normalized cfms data
message("\nLoading CFMS data w/ normalized PSMs.")
cfms_nd <- data.table::fread("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/sparklines/input_data/all_euk_cfms.concat.norm.tidy",
                             sep = ",")

# raw cfms data
message("\nLoading CFMS data w/ raw PSMs.")
cfms_rd <- data.table::fread("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/sparklines/input_data/all_euk_cfms.concat.raw.tidy",
                             sep = ",")
  
# kog <-> human <-> arath gene names
message("\nLoading euNOG annotations.")
annots <- read_tsv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/sparklines/input_data/leca_euNOGs_annotated.fmt.tsv")

# load in formal species names
message("\nLoading species information.")
species_labels <- read_csv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/sparklines/input_data/leca_experiments.csv") %>%
  janitor::clean_names() %>%
  select(code, species, species_simplified) %>% 
  mutate(code = tolower(code)) %>%
  mutate(species_simplified = str_replace(species_simplified,
                                          "Ca. Methanoperedens nitroreducens",
                                          "M. nitroreducens")) %>%
  mutate(species_short = str_remove(species_simplified, " \\(.*")) %>%
  unique()

# ordered codes based on phylogeny
species_order <- read_delim('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/euk_codes_ordered_phylo.txt', delim='\n', col_names=FALSE) %>% pull(X1)
clade_order <- c('Amorphea','Excavate','TSAR','Archaeplastida')

# -----------------------------------------------------
#################### functions ####################
# -----------------------------------------------------

get_cmplx <- function(cmplx_file, sep = ',',
                      best_exp = TRUE, sample_clades = FALSE){
  
  # read in desired complx
  cmplx <- read_delim(cmplx_file, delim = sep)
  
  cmplx_name <- cmplx %>%
    pull(granulated_cmplx_name) %>%
    unique()
  
  cmplx_ids <- cmplx %>%
    pull(ID) %>%
    unique()
  
  message("Pulling sparklines for: ", cmplx_name)
  
  cmplx_annots <- annots %>%  # get annotations
    filter(ID %in% cmplx_ids)
  
  message("IDs corresponding to this complex: \n", paste(cmplx_annots$ID, "\n"))
  
  message("Annotations:")
  print(cmplx_annots)
  
  # get representative experiment for each protein in each species if best_exp = "TRUE" (default)
  message("Extracting the complex from CFMS data...")
  if(best_exp == TRUE){
    
    message("Extracting best-sampled experiment ...")
    cmplx_rd <- cfms_rd %>%  # get complex from *ALL RAW* data
      filter(orthogroup %in% cmplx_ids) # (takes awhile)
    
    # extract best sampled experiment for each species
    cmplx_rd_ranked <- cmplx_rd %>%
      group_by(species, experiment) %>%
      tally(PSMs) %>%
      slice_max(order_by = n, n = 1) %>% 
      mutate(filter_col = paste(species, experiment, sep = " "))
      
    # extract from normalized data
    cmplx_eluts <- cfms_nd %>%  # get complex from *ALL NORMALIZED* data
      filter(orthogroup %in% cmplx_ids) %>%  # (takes awhile)
      mutate(filter_col = paste(species, experiment, sep = " ")) %>%
      filter(filter_col %in% cmplx_rd_ranked$filter_col) %>%
      mutate(line_size = 1.25)
    
  } else {
    
    cmplx_eluts <- cfms_nd %>%  # get complex from *ALL NORMALIZED* data
      filter(orthogroup %in% cmplx_ids) %>%  # (takes awhile)
      mutate(line_size = 1)
    
  }
  
  if(sample_clades == TRUE){
    
    message("Extracting subset of species from each clade ...")
    
    species_subset <- c("dicdi", "nemve", "brart", "human",  # amorphea
                        "euggr", "tryb2",  # excavate
                        "phatc", "tetts", "plakh",  # TSAR
                        "chlre", "selml", "maize", "arath")  # archaeplastida
    
    cmplx_eluts <- cmplx_eluts %>%
      filter(species %in% species_subset) %>%  ## NOTE::FOR SPECIES SUBSETTING
      mutate(line_size = 1.5)
  }
  
  message("CFMS data:")
  print(cmplx_eluts)
  return(cmplx_eluts)
  
}

fmt_df <- function(df){
  
  message("Formatting complex data frame ...")
  
  cmplx_ids <- df %>%
    pull(orthogroup) %>%
    unique()
  
  cmplx_annots <- annots %>%  # get annotations
    filter(ID %in% cmplx_ids)
  
  cmplx_fmt <- df %>%  # format for plot
    left_join(cmplx_annots, by = c("orthogroup" = "ID")) %>%
    left_join(species_labels, by = c("species" = "code")) %>%
    mutate(set = ifelse(set == "Viridiplantae", "Archaeplastida", set)) %>% 
    mutate(set = as.factor(set)) %>%
    group_by(set) %>%  # omg ->
    arrange(desc(species_short)) %>%  # this works ->
    mutate(species_short = as.factor(species_short))  # for reordering :D
  
  message("Formatted data frame:")
  print(cmplx_fmt)
  
  return(cmplx_fmt)
  
}

generate_plot <- function(df, outfile){
  
  message("Generating plot ...")
  
  theme_set(cowplot::theme_cowplot())
  pal_npg <- c("#E64B35", "#4DBBD5", "#3C5488", "#00A087",
               "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
               "#7E6148", "#B09C85")
  
  # get line size
  line_var <- df %>%
    pull(line_size) %>%
    unique()
  
  # lock in var orders
  df$set <- factor(df$set, levels = clade_order)
  df$species <- factor(df$species, levels = species_order)
  
  # generate plot
  final_plot <- ggplot(df, aes(x = fraction_id, y = PSMs)) +
    geom_line(aes(group = human_genes_fmt, color = set),
              size = line_var) +
    facet_grid(human_genes_fmt ~ species,
               switch = "y",
               scales = "free") +
    scale_color_manual(values = pal_npg) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(angle = 0)) +
    theme(strip.background = element_rect(
      color="#454545", fill="#454545", size=1.5, linetype="solid"),
      panel.spacing.x = unit(0.2, "lines"),
      strip.text.x = element_text(size = 5.5, color = "white"),
      strip.text.y = element_text(size = 12, color = "white"),
      legend.title=element_blank(),
      legend.position = "top",
      #legend.key.height= unit(5, 'cm'),
      legend.key.width= unit(2, 'cm'),
      legend.text = element_text(size=18)) +
    guides(linetype = guide_legend(override.aes = list(size = 10)))
  final_plot
  
  return(final_plot)
  
}

# -----------------------------------------------------
#################### wrapper ####################
# -----------------------------------------------------
plot_sparklines <- function(cmplx_file, sep = ',', 
                            best_exp, sample_clades,
                            outfile = NULL){
  
  cmplx <- get_cmplx(cmplx_file, sep,
                     best_exp, sample_clades)
  cmplx_fmt <- fmt_df(cmplx)
  plot <- generate_plot(cmplx_fmt)
  
  if(!is.null(outfile)){
    outfile_name <- outfile
  } else {
    outfile_name <- tools::file_path_sans_ext(cmplx_file)
    outfile_name <- paste0(outfile_name, "_sparklines")
    if(best_exp == TRUE){
      outfile_name <- paste0(outfile_name, ".sample_exp")
    }
    if(sample_clades == TRUE){
      outfile_name <- paste0(outfile_name, ".sample_clades")
    }
  }
  
  message("Saving plots to:\n", 
          paste0(outfile_name, ".pdf\n"), 
          paste0(outfile_name, ".png"))
  
  plot %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                  width = 16, height = 9, units = "in")
  plot %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                  width = 16, height = 9, units = "in")
  
  message("Done!")
  
}

plot_sparklines(cmplx_file = args$cmplx_file, outfile = args$outfile, sep = args$sep, best_exp = args$best_exp, sample_clades = args$sample_clades)
