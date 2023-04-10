library(tidyverse)
workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/ppi_predict/"
setwd(workdir)
theme_set(cowplot::theme_cowplot())
pal <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
         "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
         "#7E6148", "#B09C85")

files <- dir("feature_sweep", recursive = T, pattern = "precision_recall_*", full.names = T)

data <- data_frame(filename = files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest() %>%
  mutate(feature_set = str_extract(filename, '(?<=/).*(?=/)')) %>%
  mutate(model = str_extract(filename, '(?<=recall_).*(?=.csv)'))
data

fct_order <- c("5", "10", "25", "50", "100", "all")
data$feature_set <- factor(data$feature_set, levels = fct_order)

data %>%
  group_by(model, feature_set) %>%
  ggplot(aes(x = recall, y = precision,
             color = feature_set, group = feature_set)) +
  geom_line(size = 1.25) +
  #scale_color_manual(values = pal, name = "# features") +
  scale_color_brewer(type = "seq", palette = "YlGnBu", name = "# features") +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow=1, byrow=TRUE)) +
  geom_hline(yintercept = 0.9, linetype='dotted') +
  annotate("text", x = 0.9, y = 0.9, label = "10% FDR", vjust = -0.5) +
  facet_wrap(~model, ncol = 1)

ggsave("../../figures/pr_curve_n-feat_x_model_comparison.png", device = "png", 
       width = 6, height = 14, units = "in")


ppi_files <- dir("feature_sweep", recursive = T, pattern = "scored_interactions_fdr10*", full.names = T)

ppi_data <- data_frame(filename = ppi_files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest() %>%
  mutate(feature_set = str_extract(filename, '(?<=/).*(?=/)')) %>%
  mutate(model = str_extract(filename, '(?<=fdr10_).*(?=.csv)'))
ppi_data

options(scipen = 999)
ppi_data$feature_set <- factor(ppi_data$feature_set, levels = fct_order)
ppi_data %>%
  group_by(model, feature_set) %>%
  ggplot(aes(x = feature_set, fill = feature_set)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(type = "seq", palette = "YlGnBu", name = "# features") +
  facet_wrap(~model, ncol = 1) +
  labs(x = "", y = "# PPIs with FDR <= 10%")
  
