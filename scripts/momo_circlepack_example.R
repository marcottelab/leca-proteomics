# example input file
momo_example <- read_csv('/stor/work/Marcotte/project/Momo/rbc_figures/circlepack/category_annotated_modified_justrbc.csv') 

#label each circle with gene names
momo_label_ex <- read_csv("/stor/work/Marcotte/project/Momo/diffusion_clustering/human_proteome_uiprotID_genename.csv")
momo_label_ex <- momo_label_ex %>% 
  group_by(ID) %>% 
  mutate(num=n()) %>% 
  ungroup() %>% 
  filter(num==1)

# make border annotations
momo_bordering_ex  <- read_csv('/stor/work/Marcotte/project/Momo/rbc_figures/circlepack/category_annotated_modified_justrbc_bordering_oxidoreductase.csv')

# pick clst cols
sel <- names(momo_example)[c(3,5,8)] 

clusters_uniqued_long <- format_clusters(momo_example, sel)

graph <- make_graph_clust(clusters_uniqued_long, sel)

graph <- graph %>% activate("nodes") %>%
  left_join(momo_label_ex, by = c("name" = "ID"))%>%
  left_join(momo_bordering_ex, by = c("name"  = "ID")) 
circlepack_fxn(graph, 68)

# save output
ggsave("/stor/work/Marcotte/project/Momo/rbc_figures/circlepack/circlepackplot_category_cut358_rbconly_4cats_cb_friendly.pdf", device = "pdf", width = 6, height = 6, units = "in")