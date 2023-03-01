library(tidyverse)
theme_set(cowplot::theme_cowplot())

featselect <- read_csv("/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/feat_selection/gold.norm.102021.featselect")

models <- featselect %>%
  pull(model) %>%
  unique()

rank <- featselect %>%
  group_by(feature) %>%
  summarize(median_rank = median(rank)) %>%
  write_csv("/stor/work/Marcotte/project/rmcox/features_median_rank.allmodels.csv")


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
