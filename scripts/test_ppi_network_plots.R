library(tidyverse)
library(network)
library(igraph)
library(tidygraph)
library(ggraph)
theme_set(theme_bw(base_size = 12))
palette_pretty <- c("#0072B2","#E69F00","#009E24","#FF0000", "#5530AA")
## the custom function using Color Brewer
cols_f <- colorRampPalette(RColorBrewer::brewer.pal(8, 'Spectral'))
spectral <- as.vector(cols_f(8))

workdir <- "/stor/work/Marcotte/project/rmcox/"

# scores & complexes
ppi_scores <- read_tsv(file.path(workdir, "LECA_archive/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/scored_interactions_top30k_fdr15"))
clustering <- read_csv(file.path(workdir, "LECA_archive/ms/cfms2/cfmsflow_012022/results_lj.noribo/clustering_fdr15_noribo.fmt.csv"))

# test train status
testfile <- "/stor/work/Marcotte/project/rmcox/LECA_archive/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.test_ppis.ordered"
trainfile <- "/stor/work/Marcotte/project/rmcox/LECA_archive/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered"

test <- readr::read_csv(testfile, col_names = F) %>%
  mutate(set = "test") %>%
  rename(ID = X1)
train <- readr::read_csv(trainfile, col_names = F) %>%
  mutate(set = "train") %>%
  rename(ID = X1)

test_train <- rbind(test, train)

# combine with cfms score
ppi_scores <- left_join(ppi_scores, test_train)
ppi_scores[is.na(ppi_scores)] <- "predicted"

# ///////////////////////////////////////////////////////////////////////////////////////////

# known complexes

# filter for desired complex
complex = "Proteasome"
cluster_cut_level = 3  # between 1-9; lower = bigger clusters

cut_col <- colnames(clustering[(cluster_cut_level + 1)])

clst_filt <- clustering %>%
  select(ID, (cluster_cut_level + 1), animals, plants, tsar, excavate,
         human_gene_names_primary, human_protein_names) %>% 
  filter(grepl(complex, human_protein_names, ignore.case = TRUE)) %>%
  mutate(gene_name_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(gene_name_fmt = str_replace(gene_name_fmt, ", ", "\n"))

# get scores for desired complex
ppi_scores_filt <- ppi_scores %>%
  filter(grepl(paste(clst_filt$ID, collapse = "|"), ID, ignore.case = TRUE)) %>% 
  separate(ID, into = c("ID1", "ID2"), sep = " ")

# made node list
nodes <- clst_filt %>%
  select(ID) %>%
  rename(label = ID) %>%
  rowid_to_column("id")

# make edge list
edges <- ppi_scores_filt %>%
  rename(weight = P_1) %>% 
  left_join(nodes, by = c("ID1" = "label")) %>%
  rename(from = id) %>%
  left_join(nodes, by = c("ID2" = "label")) %>%
  rename(to = id) %>%
  select(from, to, weight, set) %>%
  na.omit()

# make (ggraph-specific) vertices list
vertices <- clst_filt %>%
  rowid_to_column("idx") %>%
  filter(idx %in% c(edges$to, edges$from))


# --------- using igraph package ---------
# library(igraph)
# ppi_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# 
# # plot network
# plot(ppi_igraph, edge.size = 0.1, layout = layout_with_graphopt)
# --------- --------- --------- -----------

# --------- using tidygraph package ---------
# reload libraries in case functions are masked by other network libraries
library(tidygraph)
library(ggraph)

ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
# or like this
# ppi_tidygraph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
# or from an igraph object
# ppi_tidygraph <- as_tbl_graph(ppi_igraph)

# the nodes tibble is "activated" by default for manipulation
# can activate edges tibble like so:
ppi_tidygraph %>%
  activate(edges) %>%
  arrange(desc(weight))

# plot network
set.seed(13)

# auto-pick a layout
ppi_network <- ggraph(ppi_tidygraph, layout = "nicely") +
  geom_edge_link(aes(width = weight),
                 alpha = 0.25) +
  geom_node_point(size = 20, aes(color = as.factor(get(cut_col))),
                  alpha = 0.9) +
  scale_edge_width_continuous(range = c(0.1, 6)) +
  geom_node_text(aes(label = gene_name_fmt),
                 repel = FALSE,
                 size = 4,
                 color = "black",
                 fontface = "bold") +
  labs(edge_width = "Interaction Score",
       color = "Cluster",
       title = complex) +
  scale_color_manual(values = palette_pretty) +
  theme_graph() +
  theme(legend.position="none") #+
  #coord_fixed()
ppi_network

ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".pdf"), ., device = "pdf", 
                      width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".png"), ., device = "png", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".tiff"), ., device = "tiff", 
                       width = 13, height = 7.5, units = "in")

clst_filt %>% write_csv(paste0(workdir, "figures/", complex, "ProtAnnots", Sys.Date(), ".csv"))

# # tree layout
# ggraph(ppi_tidygraph, 'circlepack') + 
#   geom_edge_link() + 
#   geom_node_point(aes(colour = get(cut_col)), size = 3) +
#   theme_graph() +
#   geom_node_text(aes(label = human_gene_names_primary, filter = leaf),
#                  size = 2)
# 
# # circle graph
# ggraph(ppi_tidygraph, layout = 'circlepack') + 
#   geom_node_circle(aes(fill = depth)) +
#   theme_void() +
#   theme(legend.position = "FALSE") +
#   scale_fill_viridis() +
#   geom_node_text(aes(label = human_gene_names_primary, filter = leaf), 
#                  size = 2, repel = TRUE)

# --------- --------- --------- -----------

# ///////////////////////////////////////////////////////////////////////////////////////////

# complexes w/ novel components
cut_level = 7
cmplx_group = 20
cmplx_name = "Exosome"
new_mems = c("SLC4A1AP", "ANKZF1")
  
# filter for desired complex
cut_col <- colnames(clustering[(cut_level + 1)])

clst_filt_novel <- clustering %>%
  select(ID, (cut_level + 1), animals, plants, tsar, excavate,
         human_gene_names_primary, human_protein_names) %>% 
  filter(get(cut_col) == cmplx_group) %>%
  mutate(gene_name_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(gene_name_fmt = str_replace(gene_name_fmt, ", ", "\n")) %>%
  mutate(status = case_when(grepl(paste(new_mems, collapse = "|"), human_gene_names_primary, ignore.case = TRUE) ~ "Novel",
                            TRUE ~ "Known"))

# get scores for desired complex
ppi_scores_filt <- ppi_scores %>%
  filter(grepl(paste(clst_filt_novel$ID, collapse = "|"), ID, ignore.case = TRUE)) %>% 
  separate(ID, into = c("ID1", "ID2"), sep = " ")

# made node list
nodes <- clst_filt_novel %>%
  select(ID) %>%
  rename(label = ID) %>%
  rowid_to_column("id")

# make edge list
edges <- ppi_scores_filt %>%
  rename(weight = P_1) %>% 
  left_join(nodes, by = c("ID1" = "label")) %>%
  rename(from = id) %>%
  left_join(nodes, by = c("ID2" = "label")) %>%
  rename(to = id) %>%
  select(from, to, weight, set) %>%
  na.omit()

# make (ggraph-specific) vertices list
vertices <- clst_filt_novel %>%
  rowid_to_column("idx") %>%
  filter(idx %in% c(edges$to, edges$from))

ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)

# plot network
set.seed(13)

# auto-pick a layout
ppi_network <- ggraph(ppi_tidygraph, layout = "nicely") +
  geom_edge_link(aes(width = weight),
                 alpha = 0.25) +
  geom_node_point(size = 12.8, aes(color = status)) +
  scale_edge_width_continuous(range = c(0.1, 6)) +
  geom_node_text(aes(label = gene_name_fmt),
                 repel = FALSE,
                 size = 2.3,
                 color = "white",
                 fontface = "bold") +
  labs(edge_width = "Interaction Score",
       color = "Interaction status",
       title = cmplx_name) +
  scale_color_manual(values = palette_pretty) +
  theme_graph()
ppi_network

new_mems_fmt <- paste(new_mems, collapse = "-")

ppi_network %>% ggsave(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, Sys.Date(), ".pdf"), ., device = "pdf", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, Sys.Date(), ".png"), ., device = "png", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, Sys.Date(), ".tiff"), ., device = "tiff", 
                       width = 13, height = 7.5, units = "in")

clst_filt %>% write_csv(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, "ProtAnnots", Sys.Date(), ".csv"))


# ///////////////////////////////////////////////////////////////////////////////////////////

# not human but highly conserved
nog_ids <- c("KOG0692", "KOG0560", "KOG4201")

cut_level = 6
cmplx_group = 308
cmplx_name = ""

# filter for desired complex
cut_col <- colnames(clustering[(cut_level + 1)])

clst_filt_nh <- clustering %>%
  select(ID, (cut_level + 1), animals, plants, tsar, excavate,
         human_gene_names_primary, human_protein_names) %>% 
  filter(get(cut_col) == cmplx_group) %>%
  mutate(gene_name_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(gene_name_fmt = str_replace(gene_name_fmt, ", ", "\n"))

# get scores for desired complex
ppi_scores_filt <- ppi_scores %>%
  filter(grepl(paste(nog_ids, collapse = "|"), ID, ignore.case = TRUE)) %>% 
  separate(ID, into = c("ID1", "ID2"), sep = " ")

# made node list
nodes <- clst_filt_nh %>%
  select(ID) %>%
  rename(label = ID) %>%
  rowid_to_column("id")

# make edge list
edges <- ppi_scores_filt %>%
  rename(weight = P_1) %>% 
  left_join(nodes, by = c("ID1" = "label")) %>%
  rename(from = id) %>%
  left_join(nodes, by = c("ID2" = "label")) %>%
  rename(to = id) %>%
  select(from, to, weight, set) %>%
  na.omit()

# make (ggraph-specific) vertices list
vertices <- clst_filt_nh %>%
  rowid_to_column("idx") %>%
  filter(idx %in% c(edges$to, edges$from))

ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)

# plot network
set.seed(13)

# auto-pick a layout
ppi_network <- ggraph(ppi_tidygraph, layout = "nicely") +
  geom_edge_link(aes(width = weight),
                 alpha = 0.25) +
  geom_node_point(size = 12.8, aes(color = status)) +
  scale_edge_width_continuous(range = c(0.1, 6)) +
  geom_node_text(aes(label = gene_name_fmt),
                 repel = FALSE,
                 size = 2.3,
                 color = "white",
                 fontface = "bold") +
  labs(edge_width = "Interaction Score",
       color = "Interaction status",
       title = cmplx_name) +
  scale_color_manual(values = palette_pretty) +
  theme_graph()
ppi_network

new_mems_fmt <- paste(new_mems, collapse = "-")

ppi_network %>% ggsave(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, Sys.Date(), ".pdf"), ., device = "pdf", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, Sys.Date(), ".png"), ., device = "png", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, Sys.Date(), ".tiff"), ., device = "tiff", 
                       width = 13, height = 7.5, units = "in")

clst_filt %>% write_csv(paste0(workdir, "figures/", cmplx_name, "_", new_mems_fmt, "ProtAnnots", Sys.Date(), ".csv"))



# ///////////////////////////////////////////////////////////////////////////////////////////

# cartoonized complexes

# filter for desired complex
complex = "Nuclear pore complex"
cluster_cut_level = 3  # between 1-9; lower = bigger clusters

cut_col <- colnames(clustering[(cluster_cut_level + 1)])

clst_filt <- clustering %>%
  select(ID, (cluster_cut_level + 1), animals, plants, tsar, excavate,
         human_gene_names_primary, human_protein_names) %>% 
  filter(grepl(complex, human_protein_names, ignore.case = TRUE)) %>%
  mutate(gene_name_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(gene_name_fmt = str_replace(gene_name_fmt, ", ", "\n"))

# get scores for desired complex
ppi_scores_filt <- ppi_scores %>%
  filter(grepl(paste(clst_filt$ID, collapse = "|"), ID, ignore.case = TRUE)) %>% 
  separate(ID, into = c("ID1", "ID2"), sep = " ")

# made node list
nodes <- clst_filt %>%
  select(ID) %>%
  rename(label = ID) %>%
  rowid_to_column("id")

# make edge list
edges <- ppi_scores_filt %>%
  rename(weight = P_1) %>% 
  left_join(nodes, by = c("ID1" = "label")) %>%
  rename(from = id) %>%
  left_join(nodes, by = c("ID2" = "label")) %>%
  rename(to = id) %>%
  select(from, to, weight, set) %>%
  na.omit()

# make (ggraph-specific) vertices list
vertices <- clst_filt %>%
  rowid_to_column("idx") %>%
  filter(idx %in% c(edges$to, edges$from))

ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)

# plot network
set.seed(13)

# auto-pick a layout
ppi_network <- ggraph(ppi_tidygraph, layout = "nicely") +
  geom_edge_link(aes(width = weight),
                 alpha = 0.25) +
  geom_node_point(size = 3,
                  alpha = 0.9) +
  scale_edge_width_continuous(range = c(0.1, 3)) +
  # geom_node_text(aes(label = gene_name_fmt),
  #                repel = FALSE,
  #                size = 4,
  #                color = "black",
  #                fontface = "bold") +
  labs(title = complex) +
  scale_color_manual(values = palette_pretty) +
  theme_graph() +
  theme(legend.position="none") +
  coord_fixed()
ppi_network

ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".pdf"), ., device = "pdf", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".png"), ., device = "png", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".tiff"), ., device = "tiff", 
                       width = 13, height = 7.5, units = "in")

clst_filt %>% write_csv(paste0(workdir, "figures/", complex, "ProtAnnots", Sys.Date(), ".csv"))

# # tree layout
# ggraph(ppi_tidygraph, 'circlepack') + 
#   geom_edge_link() + 
#   geom_node_point(aes(colour = get(cut_col)), size = 3) +
#   theme_graph() +
#   geom_node_text(aes(label = human_gene_names_primary, filter = leaf),
#                  size = 2)
# 
# # circle graph
# ggraph(ppi_tidygraph, layout = 'circlepack') + 
#   geom_node_circle(aes(fill = depth)) +
#   theme_void() +
#   theme(legend.position = "FALSE") +
#   scale_fill_viridis() +
#   geom_node_text(aes(label = human_gene_names_primary, filter = leaf), 
#                  size = 2, repel = TRUE)

# --------- --------- --------- -----------


# ///////////////////////////////////////////////////////////////////////////////////////////

# obligate pairs

clustering <- read_csv("/stor/work/Marcotte/project/rmcox/LECA_archive/results/leca_clustering_strongpairs.csv") %>%
  mutate(human_gene_names_primary = coalesce(human_gene_names_primary, ID))

phylo <- clustering %>%
  select(filtered_cmplx_name, animals, plants, tsar, excavate) %>%
  group_by(filtered_cmplx_name) %>%
  summarise(across(everything(), list(sum)),
            group_len = n()) %>% 
  pivot_longer(animals_1:excavate_1, 
               names_to = 'clade', 
               values_to = 'count') %>%
  mutate(conserved_interaction = as.factor(ifelse(count == group_len, 1, 0))) %>%
  mutate(clade_proper = case_when(clade == 'animals_1' ~ 'Amorphea',
                                  clade == 'plants_1' ~ 'Viridiplantae',
                                  clade == 'tsar_1' ~ 'TSAR',
                                  clade == 'excavate_1' ~ 'Excavate')) %>%
  mutate(axis = 1) %>%
  arrange(desc(filtered_cmplx_name))

ggplot(phylo, aes(x = clade_proper, y = filtered_cmplx_name)) +
  geom_dotplot(binwidth = 1, stroke = 1, binaxis = "y",
               stackdir = "center", aes(fill = conserved_interaction)) +
  scale_fill_manual(values = c("white", "#1E1E1E")) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none")

cmplx_names <- pull(clustering, filtered_cmplx_name) %>% unique()

theme_graph(base_family="sans")
  
complex <- c
cluster_cut_level <- 7  # between 1-9; lower = bigger clusters

cut_col <- colnames(clustering[(cluster_cut_level + 1)])

print(complex)
Sys.sleep(2)

clst_filt <- clustering %>%
  select(ID, (cluster_cut_level + 1), animals, plants, tsar, excavate, 
         filtered_cmplx_name, human_gene_names_primary, human_protein_names) %>% 
  #filter(filtered_cmplx_name == complex) %>%
  mutate(gene_name_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(gene_name_fmt = str_replace(gene_name_fmt, ", ", "\n")) %>% 
  mutate(gene_name_fmt = str_replace(gene_name_fmt, "NA, ", ""))

print(clst_filt)
Sys.sleep(2)

# get scores for desired complex
ppi_scores_filt <- ppi_scores %>%
  filter(grepl(paste(clst_filt$ID, collapse = "|"), ID, ignore.case = TRUE)) %>% 
  separate(ID, into = c("ID1", "ID2"), sep = " ")

# made node list
nodes <- clst_filt %>%
  select(ID) %>%
  rename(label = ID) %>%
  rowid_to_column("id")

# make edge list
edges <- ppi_scores_filt %>%
  rename(weight = P_1) %>% 
  left_join(nodes, by = c("ID1" = "label")) %>%
  rename(from = id) %>%
  left_join(nodes, by = c("ID2" = "label")) %>%
  rename(to = id) %>%
  select(from, to, weight, set) %>%
  na.omit()

# make (ggraph-specific) vertices list
vertices <- clst_filt %>%
  rowid_to_column("idx") %>%
  filter(idx %in% c(edges$to, edges$from))

# reload libraries in case functions are masked by other network libraries
library(tidygraph)
library(ggraph)

ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
# or like this
# ppi_tidygraph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
# or from an igraph object
# ppi_tidygraph <- as_tbl_graph(ppi_igraph)

# the nodes tibble is "activated" by default for manipulation
# can activate edges tibble like so:
# ppi_tidygraph %>%
#  activate(edges) %>%
#  arrange(desc(weight))

# plot network
set.seed(13)
random_fill = sample(palette_pretty, 1)
pretty_multiplied = replicate(10, palette_pretty)

# auto-pick a layout
ppi_network <- ggraph(ppi_tidygraph, layout = "nicely") +
  geom_edge_link(aes(width = weight),
                 alpha = 0.25) +
  geom_node_point(size = 20, aes(color = as.factor(get(cut_col))),
                  alpha = 0.9) +
  scale_edge_width_continuous(range = c(0.1, 6)) +
  geom_node_text(aes(label = gene_name_fmt),
                 repel = FALSE,
                 size = 4,
                 color = "black",
                 fontface = "bold") +
  labs(edge_width = "Interaction Score",
       color = "Cluster") +
  scale_color_manual(values = pretty_multiplied) +
  theme_graph() +
  theme(legend.position="none") #+
#coord_fixed()
ppi_network

ppi_network %>% ggsave(paste0(workdir, "leca_pairs_ppi_networks/", "pairs-network_", Sys.Date(), ".pdf"), ., device = "pdf", 
                       width = 10, height = 8, units = "in")
ppi_network %>% ggsave(paste0(workdir, "leca_pairs_ppi_networks/", "pairs-network_", Sys.Date(), ".png"), ., device = "png", 
                       width = 10, height = 8, units = "in")

clst_filt %>% write_csv(paste0(workdir, "leca_pairs_ppi_networks/", "pairs-network_ProtAnnots", Sys.Date(), ".csv"))

# get clade info
phylo_filt <- phylo %>%
  filter(filtered_cmplx_name == complex)

# plot clade info
theme_set(theme_minimal())
clade_dots <- ggplot(phylo, aes(x = clade_proper, y = filtered_cmplx_name)) +
  geom_dotplot(binwidth = 1, stroke = 1, binaxis = "y",
               stackdir = "center", aes(fill = conserved_interaction)) +
  scale_fill_manual(values = c("white", "#1E1E1E")) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none")
clade_dots

clade_dots %>% ggsave(paste0(workdir, "leca_pairs_ppi_networks/", complex, "phylodist_", Sys.Date(), ".pdf"), ., device = "pdf", 
                       width = 8, height = 4, units = "in")
clade_dots %>% ggsave(paste0(workdir, "leca_pairs_ppi_networks/", complex, "phylodist_", Sys.Date(), ".png"), ., device = "png", 
                       width = 8, height = 4, units = "in")

# --------- using igraph package ---------
# library(igraph)
# ppi_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# 
# # plot network
# plot(ppi_igraph, edge.size = 0.1, layout = layout_with_graphopt)
# --------- --------- --------- -----------

# --------- using tidygraph package ---------
# reload libraries in case functions are masked by other network libraries
library(tidygraph)
library(ggraph)

ppi_tidygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
# or like this
# ppi_tidygraph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
# or from an igraph object
# ppi_tidygraph <- as_tbl_graph(ppi_igraph)

# the nodes tibble is "activated" by default for manipulation
# can activate edges tibble like so:
ppi_tidygraph %>%
  activate(edges) %>%
  arrange(desc(weight))

# plot network
set.seed(13)

# auto-pick a layout
ppi_network <- ggraph(ppi_tidygraph, layout = "nicely") +
  geom_edge_link(aes(width = weight),
                 alpha = 0.25) +
  geom_node_point(size = 20, aes(color = as.factor(get(cut_col))),
                  alpha = 0.9) +
  scale_edge_width_continuous(range = c(0.1, 6)) +
  geom_node_text(aes(label = gene_name_fmt),
                 repel = FALSE,
                 size = 4,
                 color = "black",
                 fontface = "bold") +
  labs(edge_width = "Interaction Score",
       color = "Cluster",
       title = complex) +
  scale_color_manual(values = palette_pretty) +
  theme_graph() +
  theme(legend.position="none") #+
  #coord_fixed()
ppi_network

ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".pdf"), ., device = "pdf", 
                      width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".png"), ., device = "png", 
                       width = 13, height = 7.5, units = "in")
ppi_network %>% ggsave(paste0(workdir, "figures/", complex, Sys.Date(), ".tiff"), ., device = "tiff", 
                       width = 13, height = 7.5, units = "in")

clst_filt %>% write_csv(paste0(workdir, "figures/", complex, "ProtAnnots", Sys.Date(), ".csv"))

# # tree layout
# ggraph(ppi_tidygraph, 'circlepack') + 
#   geom_edge_link() + 
#   geom_node_point(aes(colour = get(cut_col)), size = 3) +
#   theme_graph() +
#   geom_node_text(aes(label = human_gene_names_primary, filter = leaf),
#                  size = 2)
# 
# # circle graph
# ggraph(ppi_tidygraph, layout = 'circlepack') + 
#   geom_node_circle(aes(fill = depth)) +
#   theme_void() +
#   theme(legend.position = "FALSE") +
#   scale_fill_viridis() +
#   geom_node_text(aes(label = human_gene_names_primary, filter = leaf), 
#                  size = 2, repel = TRUE)

# --------- --------- --------- -----------

# ///////////////////////////////////////////////////////////////////////////////////////////

