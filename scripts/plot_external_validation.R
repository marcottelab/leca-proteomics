library(tidyverse)
library(ggrepel)

workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/"
theme_set(cowplot::theme_cowplot())
pal <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
         "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
         "#7E6148", "#B09C85")
options(scipen=10000)

data <- data.table::fread(paste0(workdir, "ppi_scores_xlms_y2h_logloss.csv")) %>%
  mutate(label = ifelse(set=='xlms', 'xlms(n=912)', 'y2h(n=4610)'))

p1 <- data %>%
  group_by(set) %>% 
  ggplot(., aes(x = ppi_score, y = logloss, color = set)) +
  #geom_line(size = 2) +
  geom_point(size = 2) +
  scale_color_manual(values = pal) +
  theme(legend.position = "bottom") +
  facet_wrap(~label)

p1 %>% ggsave("leca/ppi_ml/figures/logloss_external_ppis.png", ., device = "png",
                 width = 10, height = 6, units = "in")

# p %>% ggsave("leca/ppi_ml/figures/logloss_external_ppis.pdf", ., device = "pdf",
#                  width = 6, height = 8, units = "in")

p2 <- data %>%
  group_by(set) %>% 
  ggplot(., aes(x = ppi_score, y = loglike, color = set)) +
  geom_line(size = 2) +
  #geom_point(size = 2) +
  scale_color_manual(values = pal) +
  theme(legend.position = "bottom") +
  facet_wrap(~label)

p2 %>% ggsave("leca/ppi_ml/figures/loglike_external_ppis.png", ., device = "png",
              width = 10, height = 6, units = "in")
