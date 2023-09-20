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

fct_order <- c("5", "10", "25", "50", "100", "250", "all")
data$feature_set <- factor(data$feature_set, levels = fct_order)

data %>%
  group_by(model, feature_set) %>%
  ggplot(aes(x = recall, y = precision,
             color = feature_set, group = feature_set)) +
  geom_line(size = 1.25) +
  #scale_color_manual(values = pal, name = "# features") +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow=1, byrow=TRUE)) +
  scale_color_viridis_d(option = "F", name = "# features", guide = "none", direction = -1) +
  geom_hline(yintercept = 0.9, linetype='dotted') +
  annotate("text", x = 0.9, y = 0.9, label = "10% FDR", vjust = -0.5) +
  facet_wrap(~model, nrow = 1) -> p1
p1

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
  filter(set == "predict") %>% 
  ggplot(aes(x = feature_set, fill = feature_set)) +
  geom_bar(position = "dodge") +
  scale_fill_viridis_d(option = "F", name = "# features", guide = "none", direction = -1) +
  guides(fill = guide_legend(nrow=1, byrow=TRUE)) +
  facet_wrap(~model, nrow = 1) +
  labs(x = "", y = "# predicted PPIs\nwith FDR <= 10%") -> p2

set_counts <- ppi_data %>%
  group_by(model, feature_set, set) %>%
  tally() %>%
  arrange(desc(n))

prot_counts <- ppi_data %>%
  separate(ID, c("ID1", "ID2")) %>%
  pivot_longer(ID1:ID2, values_to = "IDs") %>%
  select(-name, -ppi_score, -set) %>% 
  unique() %>%
  group_by(model, feature_set) %>%
  tally() %>%
  arrange(desc(n))

prot_counts$feature_set <- factor(prot_counts$feature_set, levels = fct_order)
prot_counts %>%
  group_by(model, feature_set) %>%
  ggplot(aes(x = feature_set, y = n, fill = feature_set)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(option = "F", name = "# features", guide = "none", direction = -1) +
  guides(fill = guide_legend(nrow=1, byrow=TRUE)) +
  facet_wrap(~model, nrow = 1) +
  labs(x = "", y = "# unique proteins in PPIs\nwith FDR <= 10%") -> p3
p3

library(patchwork)
panel <- p1 / p2 / p3 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.box.just = "center") & 
  guides(color = "none")
panel

panel %>% ggsave("../../figures/n-feat_x_model_comparison_panel_090723.png", ., device = "png",
             width = 10, height = 10, units = "in")
panel %>% ggsave("../../figures/n-feat_x_model_comparison_panel_090723.pdf", ., device = "pdf",
             width = 10, height = 10, units = "in")
