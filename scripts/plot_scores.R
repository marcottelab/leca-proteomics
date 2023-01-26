#################################################
# THIS SCRIPT WAS DEVELOPED AND TESTED BY:
# Rachael M. Cox
# rachaelcox@utexas.edu
# Last updated: 01/18/2023
#################################################

library(argparse)
library(tidyverse)
library(patchwork)
library(scales)
parser <- ArgumentParser()
parser$add_argument("-i", "--cmplx_file", dest="cmplx_file", required=TRUE,
                    help="File containing a complex (one row per protein or orthogroup), identifier column = 'ID'")
parser$add_argument("-o", "--outfile", dest="outfile",
                    help="Outfile path/name for scores figures (not required)")
parser$add_argument("-s", "--sep", dest="sep", default=",", help="Separator for infile (default = ',')")
args <- parser$parse_args()


# -----------------------------------------------------
################# global data set up ##################
# -----------------------------------------------------

# coverage data
scores <- read_csv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/apms/clustered_og_scores_subset.csv') %>%
  mutate(ID2 = coalesce(ID2, ID1))

score_cols = names(select(scores, !c(ID1, ID2, cmplx)))

# calculate min/max for each metric
min_vals = apply(scores[score_cols],2,min,na.rm=T)
max_vals = apply(scores[score_cols],2,max,na.rm=T)

# annotations
annots <- read_csv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/leca_euNOGs_plot_annots.csv') %>%
  select(-matches("human*")) %>%
  mutate(gene_names = str_replace(gene_names, " ", "\n"))

# -----------------------------------------------------
#################### functions #######################
# -----------------------------------------------------

get_data <- function(cmplx_file){

  cmplx_df <- read_csv(cmplx_file)
  cmplx_name <- cmplx_df %>%
    pull(granulated_cmplx_name) %>%
    unique()
  
  # get cluster order for prots
  cmplx_ids_ordered <- cmplx_df %>%
    select(ID, characterization_status) %>%
    unique() %>%
    mutate(fct_lvl = as.factor(row_number()))
  
  cmplx_scores <- scores %>%
    filter(cmplx == cmplx_name) %>%
    left_join(annots, by = c('ID1' = 'ID')) %>%
    left_join(cmplx_ids_ordered, by = c('ID1' = 'ID')) %>% 
    mutate(gene_names = ifelse(str_detect(characterization_status, 'Novel'),
                                   paste0('*', gene_names),
                                   gene_names)) %>%
    rename(ID1_human_gene = gene_names,
           ID1_fct_lvl = fct_lvl) %>%
    select(-characterization_status) %>% 
    left_join(annots, by = c('ID2' = 'ID')) %>%
    left_join(cmplx_ids_ordered, by = c('ID2' = 'ID')) %>%
    mutate(gene_names = ifelse(str_detect(characterization_status, 'Novel'),
                                   paste0('*', gene_names),
                                   gene_names)) %>%
    rename(ID2_human_gene = gene_names,
           ID2_fct_lvl = fct_lvl) %>%
    select(-characterization_status) %>%
    mutate(ID1_human_gene = fct_reorder(ID1_human_gene, as.numeric(ID1_fct_lvl)),
           ID2_human_gene = fct_reorder(ID2_human_gene, as.numeric(ID2_fct_lvl))) %>% 
    select(ID1, ID1_human_gene, ID2, ID2_human_gene,
           everything())
  
  return(cmplx_scores)

}

plot_scores <- function(df, score_col) {
  
  score_min = round(min_vals[[score_col]], digits=1)
  score_max = round(max_vals[[score_col]], digits=2)
  
  maxmin = paste0('(min=', score_min,', max=', score_max,')')
  title_text = paste(score_col, maxmin, sep='\n')
  
  theme_set(theme_bw())
  pal = c('#e6e172', '#e64b35')
  p <- ggplot(df, aes(x = ID1_human_gene, y = ID2_human_gene,
                           fill = get(score_col))) +
    geom_tile(alpha = 0.85, linejoin = "round") +
    geom_text(aes(label = ifelse(is.na(get(score_col)), "", 
                                 round(get(score_col), digits=2))),
              size = 3) +
    scale_fill_gradient(low=pal[1], high=pal[2],
                        limits = c(as.integer(score_min), score_max), na.value = "white",
                        breaks = c(as.integer(score_min), score_max),
                        name = paste0(score_col,':')) +
    theme(legend.margin = margin(0, 0, 0, 0),
          #legend.key.size = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.key.height = unit(0.25, "cm"),
          legend.spacing.x = unit(0, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.position = "top",
          #legend.title = element_text(size = 4),
          #legend.text = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill="white", colour="white"),
          plot.title = element_text(size = 6),
          plot.margin = unit(c(0,0.025,0.025,0), "cm")) +
    guides(fill = guide_colorbar(title.position = "top"))
  
  
  # get cmplx size
  ids1 <- pull(df, ID1) %>% unique
  ids2 <- pull(df, ID2) %>% unique
  num_units <- length(unique(append(ids1, ids2)))
  
  if(num_units >= 15){

    p$layers[[2]] <- NULL
    p <- p +
      theme(text = element_text(size=6.75),
            legend.key.size = unit(1.5, "cm"))
    
  }
  
  return(p)
  
}

plot_panel <- function(df, score_cols){
  
  num_scores = length(score_cols)
  
  p1 <- plot_scores(df = df, score_col = score_cols[1])
  pn <- lapply(score_cols[2:num_scores], plot_scores, df = df)
  print(length(pn))
  
  # only keep y axis tick labels on first plot
  for(i in seq(1, length(pn), 1)){
    pn[[i]] <- pn[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  panel = p1 + pn + plot_layout(nrow = 1, ncol = num_scores)
  panel
  return(panel)
  
}

# -----------------------------------------------------
#################### main #######################
# -----------------------------------------------------

# test code
# cmplx_file <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/Aminoacyl-tRNA_synthetase_multienzyme_complex.csv"
# score_df <- get_data(cmplx_file = cmplx_file)

score_df <- get_data(cmplx_file = args$cmplx_file)

score_plots <- plot_panel(score_df, score_cols)
score_plots

if(!is.null(args$outfile)){
  outfile_name <- args$outfile
} else {
  outfile_name <- tools::file_path_sans_ext(args$cmplx_file)
  outfile_name <- paste0(outfile_name, "_scores")
}

message("Saving plots to:\n", 
        paste0(outfile_name, ".pdf\n"), 
        paste0(outfile_name, ".png"))  

score_plots %>% ggsave(paste0(outfile_name, ".pdf"), ., device = "pdf", 
                       width = 16, height = 4, units = "in")
score_plots %>% ggsave(paste0(outfile_name, ".png"), ., device = "png", 
                       width = 16, height = 4, units = "in")

message("Done!")
