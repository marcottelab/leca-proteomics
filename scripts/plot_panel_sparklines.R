#################################################
# THIS SCRIPT WAS DEVELOPED AND TESTED BY:
# Rachael M. Cox
# rachaelcox@utexas.edu
# Last updated: 01/18/2023
#################################################
library(argparse)
library(tidyverse)
library(scales)

# normalized cfms data
message("\nLoading CFMS data w/ normalized PSMs.")
cfms_nd <- data.table::fread("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/elutions/all_euk_cfms.concat.norm.tidy",
                             sep = ",")

# raw cfms data
message("\nLoading CFMS data w/ raw PSMs.")
cfms_rd <- data.table::fread("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/elutions/all_euk_cfms.concat.raw.tidy",
                             sep = ",")

# annotations
annots <- read_csv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/leca_euNOGs_plot_annots.csv') %>%
  mutate(gene_names = str_replace(gene_names, "\n", "/"),
         gene_names = str_replace(gene_names, " ...", ""))

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
  unique() %>%
  filter(!species %in% c("Arabidopsis thaliana mutant ttg1-1", 
                          "Plasmodium falciparum (isolate 3D7)"))

# ordered codes based on phylogeny
species_order <- read_delim('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/euk_codes_ordered_phylo.txt', delim='\n', col_names=FALSE) %>% pull(X1)
clade_order <- c('Amorphea','Excavate','TSAR','Archaeplastida')

# -----------------------------------------------------
#################### functions ####################
# -----------------------------------------------------

get_cmplx <- function(cmplx, best_exp = TRUE, n_best = 1, 
                      sample_clades = FALSE, species_list = NULL){
  
  # get cluster order for prots
  cmplx_ids_ordered <- cmplx %>%
    select(ID) %>%
    unique() %>%
    mutate(fct_lvl = row_number())
  
  # get annotations
  cmplx_annots <- annots %>%
    filter(ID %in% cmplx_ids_ordered$ID)
  
  message("IDs corresponding to this complex: \n", paste(cmplx_annots$ID, "\n"))
  
  message("Annotations:")
  print(cmplx_annots)
  
  # get representative experiment for each protein in each species if best_exp = "TRUE" (default)
  message("Extracting the complex from CFMS data...")
  if(best_exp == TRUE){
    
    message("Extracting best-sampled experiment ...")
    cmplx_rd <- cfms_rd %>%  # get complex from *ALL RAW* data
      filter(orthogroup %in% cmplx_ids_ordered$ID) # (takes awhile)
    
    # extract best sampled experiment for each species
    cmplx_rd_ranked <- cmplx_rd %>%
      group_by(species, experiment) %>%
      tally(PSMs) %>%
      slice_max(order_by = n, n = n_best) %>% 
      mutate(filter_col = paste(species, experiment, sep = " "))
    
    # extract from normalized data
    cmplx_eluts <- cfms_nd %>%  # get complex from *ALL NORMALIZED* data
      filter(orthogroup %in% cmplx_ids_ordered$ID) %>%  # (takes awhile)
      mutate(filter_col = paste(species, experiment, sep = " ")) %>%
      filter(filter_col %in% cmplx_rd_ranked$filter_col) %>%
      mutate(line_size = 1)
    
  } else {
    
    cmplx_eluts <- cfms_nd %>%  # get complex from *ALL NORMALIZED* data
      filter(orthogroup %in% cmplx_ids_ordered$ID) %>%  # (takes awhile)
      mutate(line_size = 0.8)
    
  }
  
  if(sample_clades == TRUE){
    
    message("Extracting subset of species from each clade ...")
    
    if(!is.null(species_list)){
      
      species_subset <- unlist(strsplit(species_list, ","))
      print(species_subset)
      
      
    } else {
      
      species_subset <- c("nemve", "brart", "strpu", "human",  # amorphea
                          "euggr",  "tryb2", # excavate
                          "phatc", "tetts", "plakh",  # TSAR
                          "chlre", "selml", "maize", "arath")  # archaeplastida
      
      species_subset <- c("nemve", "brart", "human",  # amorphea
                          "euggr",  # excavate
                          "phatc", "tetts",  # TSAR
                          "chlre", "cocnu", "arath")  # archaeplastida
      
    }
    
    cmplx_eluts <- cmplx_eluts %>%
      filter(species %in% species_subset) %>%  ## NOTE::FOR SPECIES SUBSETTING
      mutate(line_size = 1.5)
  }
  
  # retain cluster order for subunits
  cmplx_eluts_out <- cmplx_eluts %>%
    unique() %>% 
    left_join(cmplx_ids_ordered, by = c("orthogroup" = "ID"))
  
  message("CFMS data:")
  print(cmplx_eluts_out)
  return(cmplx_eluts_out)
  
}

fmt_df <- function(cmplx_sparklines, annotation_df, species_df){
  
  message("Formatting complex data frame ...")
  
  sparks_annotated <- left_join(cmplx_sparklines, annotation_df, by=c("orthogroup"="ID"))
  
  plot_df <- sparks_annotated %>%  # format for plot
    left_join(species_labels, by = c("species" = "code")) %>%
    mutate(set = ifelse(set == "Viridiplantae", "Archaeplastida", set)) %>%
    # mutate(gene_names_nov = ifelse(str_detect(final_status, 'Novel'),
    #                                paste0('*', gene_names),
    #                                gene_names)) %>% 
    mutate(set = as.factor(set)) %>% 
    mutate(gene_names = fct_reorder(gene_names, fct_lvl)) # maintain prot cluster order
    # mutate(gene_names_nov = fct_reorder(gene_names_nov, fct_lvl))
  
  message("Formatted data frame:")
  print(plot_df)
  
  return(plot_df)
  
}

generate_plot <- function(plot_df, outfile){
  
  message("Generating plot ...")
  
  theme_set(cowplot::theme_cowplot())
  pal_npg <- c("#E64B35", "#4DBBD5", "#3C5488", "#00A087",
               "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
               "#7E6148", "#B09C85")
  
  # get line size
  line_var <- plot_df %>%
    pull(line_size) %>%
    unique()
  
  # lock in var orders
  plot_df <- plot_df %>% mutate(species = toupper(species))
  plot_df$set <- factor(plot_df$set, levels = clade_order) # major clade order
  plot_df$species <- factor(plot_df$species, levels = toupper(species_order)) # phylogenetic order
  
  print(plot_df)
  # generate plot
  final_plot <- plot_df %>%
    ggplot(aes(x = fraction_id, y = PSMs)) +
    geom_line(aes(group = gene_names, color = set),
              size = 1) +
    facet_grid(gene_names ~ species,
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
      panel.spacing.x = unit(0.15, "lines"),
      strip.text.x = element_text(size = 11, color = "white", face="bold"),
      strip.text.y = element_text(size = 10, color = "white", face="bold"),
      legend.title=element_blank(),
      legend.position = "top",
      #legend.key.height= unit(5, 'cm'),
      legend.key.width= unit(2, 'cm'),
      legend.text = element_text(size=18),
      plot.margin = unit(c(0,0.5,0.5,0), "cm")) +
    guides(linetype = guide_legend(override.aes = list(size = 10)))
  final_plot
  
  # get cmplx size
  num_units <- length(pull(plot_df, orthogroup) %>% unique)
  if(num_units <= 5){
    
    final_plot <- final_plot +
      theme(strip.text.x = element_text(size = 9, color = "white"),
            strip.text.y = element_text(size = 18, color = "white"))
  }
  
  return(final_plot)
  
}

# -----------------------------------------------------
#################### wrapper ####################
# -----------------------------------------------------
plot_sparklines <- function(cmplx_file, sep = ',', 
                            best_exp, sample_clades, species_list,
                            outfile = NULL){
  
  cmplx <- get_cmplx(cmplx_file, sep, best_exp, sample_clades, species_list)
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
  
  # get cmplx size
  num_units <- length(pull(cmplx, orthogroup) %>% unique)
  
  # save plots based on complex size
  if(num_units <= 5){
    
    print(paste0("# subunits = ", num_units))
    message("Adjusting plot dimensions ...")
    hval = 1.5+num_units
    
    plot %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                    width = 16, height = hval, units = "in")
    plot %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                    width = 16, height = hval, units = "in")
    
  } else {
    
    plot %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                    width = 16, height = 9, units = "in")
    plot %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                    width = 16, height = 9, units = "in")
    
    message("Done!")
    
  }
  
}
# 
# plot_sparklines(cmplx_file = args$cmplx_file, outfile = args$outfile, sep = args$sep, 
#                 best_exp = args$best_exp, sample_clades = args$sample_clades, 
#                 species_list = args$species_list)


df <- read_csv("ppi_ml/data/panel_sparklines/panel_sparklines.csv") %>%
  filter(group_num==0)
cmplx_sparks <- get_cmplx(df, best_exp=TRUE, n_best=2, sample_clades=TRUE)
cmplx_annots <- fmt_df(cmplx_sparks, annots)

plot <- generate_plot(cmplx_annots) + theme(legend.position = "none")
plot

plot %>% ggsave("ppi_ml/figures/panel_sparklines.pdf", ., device = "pdf", 
                width = 9, height = 12, units = "in")
plot %>% ggsave("ppi_ml/figures/panel_sparklines.png", ., device = "png", 
                width = 9, height = 12, units = "in")
