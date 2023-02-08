library(tidyverse)

workdir <- "/stor/work/Marcotte/project/rmcox/leca/"
setwd(workdir)

theme_set(cowplot::theme_cowplot())
pal <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
         "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
         "#7E6148", "#B09C85")

pr1f <- "ppi_ml/results/pr_curves/pr_cfms_only_rfe"
pr2f <- "ppi_ml/results/pr_curves/pr_cfms_plus_apms_rfe"
pr3f <- "ppi_ml/results/pr_curves/pr_cfms_apms_fixed_rfe"


pr1 <- readr::read_csv(pr1f) %>%
  #filter(label == "test") %>%
  mutate(set = "CFMS only") %>%
  #select(Recall, Precision, set)
  select(Recall, Precision, set, label)

pr2 <- readr::read_csv(pr2f) %>%
  #filter(label == "test") %>%
  mutate(set = "CFMS + APMS") %>%
  #select(Recall, Precision, set)
  select(Recall, Precision, set, label)

pr3 <- readr::read_csv(pr3f) %>%
  #filter(label == "test") %>%
  mutate(set = "CFMS + APMS (fixed; +574123 PPI scores)") %>%
  #select(Recall, Precision, set)
  select(Recall, Precision, set, label)

all_pr <- rbind(pr1, pr2, pr3)

ggplot(all_pr, aes(x = Recall, y = Precision, color = label, group = label)) +
  geom_line(size = 2) +
  scale_color_manual(values = pal, name = "") +
  theme(legend.position = "top") +
  facet_wrap(~set, ncol=1)
