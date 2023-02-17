#################################################
# THIS SCRIPT WAS DEVELOPED AND TESTED BY:
# Rachael M. Cox
# rachaelcox@utexas.edu
# Last updated: 02/16/2023
#################################################

library(argparse, warn.conflicts=FALSE)
library(tidyverse, warn.conflicts=FALSE)
library(scales, warn.conflicts=FALSE)
library(igraph, warn.conflicts=FALSE)
library(tidygraph, warn.conflicts=FALSE)
library(ggraph, warn.conflicts=FALSE)

set.seed(13)

parser <- ArgumentParser()
parser$add_argument("-i", "--cmplx", dest="cmplx", required=TRUE,
                    help="CSV file containing a single unique complex (one row per protein or orthogroup), identifier column = 'ID'")
parser$add_argument("-s", "--scores", dest="scores", required=TRUE,
                    help="CSV file containing pairwise probability scores (e.g., output from a ML pipeline)")
parser$add_argument("-a", "--annotations", dest="annotations", required=TRUE,
                    help="CSV file containing protein annotations; protein ID column must match the complex file ID column")
parser$add_argument("-o", "--outfile", dest="outfile", required=TRUE,
                    help="Outfile path/name for PPI network figures")
parser$add_argument("--scan_layouts", action="store_true", default=FALSE, help="Whether to scan multiple graph layouts")
args <- parser$parse_args()

# ------------------------------
# functions
# ------------------------------

fmt_cmplx_data <- function(cmplx_file, annot_file){
  
  # load annotations
  annots <- read_csv(annot_file)
  
  # load cmplx data
  df_cmplx <- read_csv(cmplx_file)
  
  # format data
  cmplx_fmt <- df_cmplx %>%
    select(ID, matches("cut*"), granulated_cmplx_name, characterization_status,
           animals, excavate, tsar, plants) %>%
    mutate(n_clades = rowSums(across(c(animals:plants)))) %>%
    left_join(annots) %>%
    mutate(status_fmt = case_when(characterization_status == "Known" ~ "Known",
                                  characterization_status == "Novel member" ~ "Uncharacterized",
                                  characterization_status == "Uncharacterized" ~ "Uncharacterized"))
  
  return(cmplx_fmt)
  
}

extract_ppi_scores <- function(scores_file, cmplx_fmt){
  
  # load pairwise score data
  df_scores <- read_csv(scores_file)
  
  # get pairwise scores for complex members
  cmplx_scores <- df_scores %>%
    mutate(ID = paste(ID1, ID2, " ")) %>% 
    filter(grepl(paste(cmplx_fmt$ID, collapse = "|"), ID, ignore.case = TRUE))
  
  print("# ppi scores extracted:")
  print(nrow(cmplx_scores))
  
  return(cmplx_scores)
  
}

make_nodes_edges <- function(cmplx_fmt, cmplx_scores){
  
  # made node list
  nodes <- cmplx_fmt %>%
    select(ID) %>%
    rename(label = ID) %>%
    rowid_to_column("id")
  
  # make edge list
  edges <- cmplx_scores %>%
    rename(weight = P_1) %>% 
    left_join(nodes, by = c("ID1" = "label")) %>%
    rename(from = id) %>%
    left_join(nodes, by = c("ID2" = "label")) %>%
    rename(to = id) %>%
    select(from, to, weight, set) %>%
    na.omit()
  
  # make (ggraph-specific) vertices list
  vertices <- cmplx_fmt %>%
    rowid_to_column("idx") %>%
    filter(idx %in% c(edges$to, edges$from))
  
  # store output in list
  # because R is dumb & can only output 1 object at a time
  graph_info <- list("nodes" = nodes, "edges" = edges, "vertices" = vertices)
  
  return(graph_info)
  
}

plot_network <- function(edges, vertices, node_size, text_size, graph_layout='circle'){
  
  # plot params
  theme_set(theme_bw(base_size = 6))
  glim=text_size/2
  ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
  lay = create_layout(ppi_tidygraph, layout = graph_layout)
  
  # ggraph
  ppi_network <- ggraph(lay) +
    geom_edge_link(aes(width = weight), alpha = 0.4, linejoin = "round") +
    geom_node_point(aes(size = as.factor(n_clades)), col = "white") +
    geom_node_point(aes(fill = as.factor(status_fmt),
                        alpha = as.factor(n_clades), size = as.factor(n_clades)),
                    color = "black", shape = 21, stroke = 2) +
    scale_edge_width_continuous(range = c(0.5, 5)) +
    geom_node_text(aes(label = gene_names),
                   repel = FALSE,
                   size = text_size,
                   color = "black",
                   fontface = "bold",
                   lineheight = 0.75) +
    labs(edge_width = "PPI Score",
         fill = "Association",
         size = "# Clades") +
    scale_fill_manual(values = c("Known" = "#0072B2", "Uncharacterized" = "#F0E442")) +
    scale_size_manual(values = c(node_size+1, node_size*2, node_size*3), breaks = c(2, 3, 4)) +
    scale_alpha_manual(values = c(0.35, 0.5, 0.75), breaks = c(2, 3, 4), guide = "none") +
    scale_edge_width(limits = c(0.2, 1)) +
    xlim(-glim, glim) +
    ylim(-glim, glim) +
    theme(legend.position="bottom", legend.margin=margin()) +
    #theme(legend.position="none") +
    theme_graph()
  
  if(graph_layout == 'circle'){
    
    ppi_network <- ppi_network +
      theme(legend.position="none")
  }
    
  
  return(ppi_network)
  
}

save_plots <- function(network_plot, num_subunits, outname, suffix=NULL){
  
  outname <- paste0(outname, "_", num_subunits, "mem")
  
  if(!is.null(suffix)){
    outname = paste0(outname, '_', suffix)
  }
  
  p <- network_plot #+
    #annotate("text", x = -2, y = -2, label = paste0("# subunits = ", num_subunits))
  
  p %>% ggsave(paste0(outname, ".png"), ., device = "png", 
                         width = 10, height = 10, units = "in")
  p %>% ggsave(paste0(outname, ".pdf"), ., device = "pdf", 
                         width = 10, height = 10, units = "in")
  
}


# ------------------------------
# main
# ------------------------------

cmplx_data <- fmt_cmplx_data(args$cmplx, args$annotations)
cmplx_ppi_scores <- extract_ppi_scores(args$scores, cmplx_data)
cmplx_graph <- make_nodes_edges(cmplx_data, cmplx_ppi_scores)

cmplx_size = nrow(cmplx_data)
if(cmplx_size >= 28){
  node_size = 6
  text_size = 2
} else if(cmplx_size >= 14){
  node_size = 8
  text_size = 3
} else {
  node_size = 15
  text_size = 4
}

if(args$scan_layouts == TRUE){
  
  scan_fxn <- function(layout){
    
    plot <- plot_network(cmplx_graph$edges, cmplx_graph$vertices, node_size, text_size, layout)
    save_plots(plot, cmplx_size, args$outfile, suffix = layout)
    
  }
  
  layout_options = c('circle', 'kk', 'stress')
  lapply(layout_options, scan_fxn)
  
} else {
  
  plot <- plot_network(cmplx_graph$edges, cmplx_graph$vertices, node_size, text_size)
  save_plots(plot, cmplx_size, args$outfile)
  
}





