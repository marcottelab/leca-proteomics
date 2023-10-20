#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (complex name)", call.=FALSE)
} else if (length(args)==1) {
  cmplx2label = args[1]
  outfile = paste0(cmplx2label, '.png')
  outfile = gsub(" ", "_", outfile)
}

library(tidyverse)
library(network)
library(igraph)
library(tidygraph)
library(ggraph)
theme_set(theme_bw(base_size = 12))

pal_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
             "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
             "#7E6148", "#B09C85")
workdir <- "/stor/work/Marcotte/project/rmcox/leca/"
cut_choice_status_col = 'final_status'

# ------------------
# functions
# ------------------

format_clst <- function(clusters, level_vect){
  # get clusters greater than one big
  highest_level <- level_vect[1]
  clusters_filt <- clusters %>%
    add_count(!!sym(highest_level)) %>%
    filter(n > 1) %>%
    select(-n)
  # make each cluster number unique by adding on the header
  clusters_uniqued_long <- clusters_filt %>%
    gather(clusterset, clusternum, -ID) %>%
    mutate(clusterid = paste0(clusternum, clusterset)) %>% 
    select(-clusternum)
  return(clusters_uniqued_long)
}

make_edge <- function(n, clu, cols){
  startcol = cols[n]
  endcol = cols[n + 1]
  clu %>% select(!!startcol, !!endcol) %>% rename(from = !!startcol, to = !!endcol)
}

graph_clst <- function(clusters_uniqued_long, level_vect) {
  # make wide
  clusters_uniqued <- clusters_uniqued_long %>%
    spread(clusterset, clusterid)
  # set origin  
  clusters_uniqued$link <- "origin"
  # get origin col, cut cols, and ID col
  clusters_uniqued_sel <- clusters_uniqued %>%
    select(link, {{ level_vect }}, ID)
  # get col names
  cols <- names(clusters_uniqued_sel)
  # generate edges
  edges <-
    map_df(seq(1, length(cols) - 1), make_edge, clusters_uniqued_sel, cols) %>%
    unique()
  # generate graph
  mygraph2 <- as_tbl_graph(edges)
  return(mygraph2)
}

format_labels <- function(clusters, annot, level_vect){
  
  # First get only selected clustering levels
  clusterid_annot <- clusters %>% left_join(annot, by = "ID") %>%
    filter(clusterset %in% {{ level_vect }}) %>%
    select(clusterset, clusterid, complex_label) %>%
    
    # Count numbers of members per cluster set
    group_by(clusterid) %>%
    mutate(unique_types = n_distinct(complex_label)) %>%
    add_count(clusterset, complex_label, sort = TRUE) %>%
    
    group_by(complex_label) %>%
    filter(unique_types ==  min(unique_types)) %>%
    mutate(label = case_when(row_number() == 1 ~ complex_label)) %>%
    filter(!is.na(label))  %>%
    group_by(clusterid) %>%
    summarize(label = paste(unique(label), collapse = ', ')) %>% 
    ungroup()
  
  return(clusterid_annot)
}

annotate_graph <- function(graph, clusterid_annot, annot){
  graph_annot <- graph %>% activate('nodes') %>%
    left_join(clusterid_annot, by = c("name"  = "clusterid")) %>%
    left_join(annot, by = c("name"  = "ID")) 
  return(graph_annot)
  
}

plot_graph <- function(graph, labels = TRUE){
  
  circlepack_plants_highlight <- ggraph(graph, layout = 'circlepack') + 
    geom_node_circle(aes(  fill = set), size = 0.1, color = "grey50") +
    scale_fill_manual(values = c(palette[4],  "yellow"), na.value = "white") + #light green "#c5f7d1
    theme(legend.key.size = unit(0.2, "lines"),
          legend.position = "bottom",
          legend.text = element_text(size = 4)) +
    theme_void() + 
    NULL
  
  if(labels == TRUE){
    # label nodes
    circlepack_plants_highlight <- circlepack_plants_highlight + 
      geom_node_text( aes( label= label), repel = TRUE,  segment.color = "grey40", size = 2, force = 0.5 ) +
      NULL
  }
  
  circlepack_plants_highlight_nolegend <- circlepack_plants_highlight + 
    theme(legend.position="FALSE")  
  return(circlepack_plants_highlight_nolegend)
}

circlepack_fxn <- function(graph_annot, seed){
  
  
  set.seed(seed = seed)
  base_layout_circlepack <- ggraph(graph_annot, layout = 'circlepack')
  
  base_layout_circlepack_format <- base_layout_circlepack +
    geom_node_circle(aes(fill = cmplx_label)) +  # node fill variable
    scale_fill_manual(values = c("#3C5488","#4DBBD5", "#E64B35", "#F0E442", "grey90", "#00A087"),
                      na.translate=F) + # node fill colors
    scale_color_manual(values = c("grey50") ,na.value = "grey50") +
    scale_size_manual(values = c(1), na.value = 0.5) + # set carefully; big nums freeze RStudio
    theme_void() +
    theme(legend.title = element_blank(),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  return(base_layout_circlepack_format)
}

# ------------------
# main
# ------------------

# clusters & complexes
clst_annot_file <- paste0(workdir, "ppi_ml/results/walktrap/LinearSVC_100feats_fdr10_4steps_nochloro_dynamic_algo_labels_092623.csv") 

clst_annot <- read_csv(clst_annot_file) %>%
  mutate(human_family_size = str_count(human_gene_names_primary, ',')+1) %>%
  mutate(human_family_size = ifelse(is.na(human_family_size), 0, human_family_size)) %>% 
  #mutate(human_genes_fmt = gsub(";", ",", human_gene_names_primary, fixed=TRUE)) %>%
  mutate(human_genes_fmt = str_replace_all(human_gene_names_primary, pattern = ";", replacement = ",")) %>% 
  mutate(human_genes_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(human_genes_fmt = str_replace(human_genes_fmt, ", ", "\n")) %>%
  mutate(genes_fmt = ifelse(human_family_size != 0 & human_family_size <= 5, human_genes_fmt, ID)) %>%
  rename(characterization_status = cut_choice_status_col) %>% 
  select(ID, human_gene_names_primary, human_genes_fmt, human_family_size, genes_fmt,
         matches("cut_*"), granulated_cmplx_name, characterization_status)

# format data
pdf <- clst_annot %>%
  filter(!is.na(characterization_status),
         characterization_status != "Unclustered") %>% 
  mutate(characterization_status_fmt = case_when(characterization_status == "Known" ~ "Known protein complex",
                                                 characterization_status == "Novel association" ~ "Novel association",
                                                 characterization_status == "Uncharacterized" ~ "Uncharacterized PPI",
                                                 TRUE ~ characterization_status)) 

legend_order = c("Known protein complex","Novel association","Uncharacterized PPI",
                 "Large heterogeneous complex")

pdf$characterization_status_fmt <- factor(pdf$characterization_status_fmt,
                                          levels = legend_order)

# identify complex that will be labeled
pdf_lab <- pdf %>%
  mutate(cmplx_label = ifelse(grepl(cmplx2label, granulated_cmplx_name), granulated_cmplx_name, NA))

# get cut cols
cuts <- names(clst_annot[grepl("cut", names(clst_annot))])
cuts_sel <- cuts[c(6,8,9,10)]

# format chosen clusters
clusters_uniqued_long <- format_clst(pdf_lab, cuts_sel)

# generate graph objct
graph <- graph_clst(clusters_uniqued_long, cuts_sel)

# join labels
graph <- graph %>% activate("nodes") %>%
  left_join(pdf_lab, by = c("name" = "ID"))

message(sprintf("Plotting circle plot w/ %s highlighted ...", cmplx2label))
# plot graph object
circlepack_fxn(graph, 13)

# save output
message(sprintf("Saving %s ...", outfile))
figure_dir <- paste0(workdir, "ppi_ml/figures/cplot_cmplx_labels/")
ggsave(paste0(figure_dir, outfile), device = "png", width = 10, height = 10, units = "in")

message("Done!")