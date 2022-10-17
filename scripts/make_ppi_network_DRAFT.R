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
  geom_node_point(size = 18, aes(color = as.factor(get(cut_col))),
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