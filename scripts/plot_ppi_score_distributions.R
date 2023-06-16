library(tidyverse)
library(ggrepel)

workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/ppi_predict/feature_sweep"
theme_set(cowplot::theme_cowplot())
pal <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
         "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
         "#7E6148", "#B09C85")

files <- dir(workdir, recursive = T, pattern = "scored_interactions_fdr.*", full.names = T)

data <- data_frame(filename = files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest(cols = c(contents)) %>%
  mutate(feature_set = str_extract(filename, '(?<=feature_sweep/).*(?=/scored)')) %>%
  mutate(model = str_extract(filename, '(?<=fdr10_).*(?=.csv)'))
data

fct_order <- c("5", "10", "25", "50", "100", "250", "all")
data$feature_set <- factor(data$feature_set, levels = fct_order)

ggplot(data, aes(x = model, y = ppi_score, fill = feature_set)) +
  geom_boxplot() +
  scale_fill_brewer(type = "seq", palette = "YlGnBu", name = "# features") +
  labs(x = "", y = "PPI Scores (10% FDR)") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow=1, byrow=TRUE))
