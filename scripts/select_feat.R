featselect <- read_csv("/project/rmcox/LECA/ms/cfms2/cfmsflow/feat_selection/gold.norm.102021.featselect")

models <- featselect %>%
  pull(model) %>%
  unique()

theme_set(cowplot::theme_cowplot())
featselect_plot <- featselect %>%
  ggplot(., aes(x=rank, y = score, label = feature)) +
  geom_vline(xintercept = 49, color = 'lightblue2', size=2) +
  geom_line() +
  geom_point() +
  xlab("Feature Rank") +
  ylab("Feature Importance (Gini)") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.01, 0.02), limits = c(0,0.02)) +
  facet_wrap(~model)
featselect_plot

#featselect_etc_plot  %>% ggsave("figures/featselect_etc_plot.png", .,  device = "png", height = 2, width = 2.5, units = "in")
#featselect_etc_plot  %>% ggsave("figures/featselect_etc_plot.pdf", .,  device = "pdf", height = 2, width = 2.5, units = "in")

# featselect %>%
#   filter(model == 'ExtraTreesClassifier') %>% 
#   filter(rank <= 49) %>%
#   select(feature) %>%
#   write_csv("/project/rmcox/LECA/ms/cfms2/cfmsflow/feat_selection/gold_top100_feats_etc.norm.txt", col_names = FALSE)
  

fts <- data.frame()
for(i in models) {
  featselect %>%
    filter(model == i) %>% 
    filter(rank <= 49) %>%
    select(feature) %>%
    write_csv(sprintf("/project/rmcox/LECA/ms/cfms2/cfmsflow/feat_selection/top50_feats_%s.gold.norm.csv", i), col_names = FALSE)
}

data_path <- "/project/rmcox/LECA/ms/cfms2/cfmsflow/feat_selection/"   # path to the data
featfiles <- dir(data_path, pattern = "top50_feats.*norm*") # get file names

feat_overlap <- files %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(file.path(data_path, .), col_names = FALSE)) %>% 
  reduce(rbind) %>%
  group_by(X1) %>%
  tally() %>%
  arrange(desc(n)) %>%
  head(50)
feat_overlap

write_csv(feat_overlap, paste0(data_path, 'top50_feats.allmodels.gold.norm.csv'))
