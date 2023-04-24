library(tidyverse)
library(patchwork)

dir <- '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cross_val'
files <- list.files(dir, pattern="precision_recall*", full.names = T)

fold_exp = '\\d(?=.csv)'
model_exp = '(?<=recall_).+(?=\\d)'

data <- data_frame(filename = files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest(cols = c(contents)) %>%
  mutate(fold = as.numeric(str_extract(filename, fold_exp))) %>%
  mutate(model = str_extract(filename, model_exp))
data

ap_files <- list.files(dir, pattern="*avg_precision.csv", full.names = T)
ap_data <- data_frame(filename = ap_files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest() %>%
  mutate(model = str_extract(filename, '(?<=cross_val/).+(?=_GroupKFold)')) %>%
  select(-filename)

pdf <- data %>%
  left_join(ap_data, by=c("fold" = "fold", "model" = "model")) %>%
  mutate(ap_label = paste0("(AP=", average_precision, ")")) %>% 
  mutate(legend_label = paste0("GroupKFold ",fold," ",ap_label))

theme_set(theme_bw())
palette_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85")

p1 <- pdf %>%
  filter(model == "ExtraTreesClassifier") %>% 
  group_by(model, fold) %>%
  ggplot(aes(x = recall, y = precision, color = legend_label)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = palette_npg, name = "") +
  theme(legend.position = c(0.2, 0.2),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  facet_wrap(~model)

p2 <- pdf %>% 
  filter(model == "SGDClassifier") %>% 
  group_by(model, fold) %>%
  ggplot(aes(x = recall, y = precision, color = legend_label)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = palette_npg, name = "") +
  theme(legend.position = c(0.3, 0.2),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  facet_wrap(~model)

p2 <- pdf %>% 
  filter(model == "LinearSVC") %>% 
  group_by(model, fold) %>%
  ggplot(aes(x = recall, y = precision, color = legend_label)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = palette_npg, name = "") +
  theme(legend.position = c(0.3, 0.2),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  facet_wrap(~model)

panel = p1 + p2 + p3
panel
