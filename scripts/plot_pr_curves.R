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
pr4f <- "ppi_ml/results/pr_curves/pr_new_pipeline"


pr1 <- readr::read_csv(pr1f) %>%
  filter(label == "test") %>%
  mutate(set = "CFMS only\n(includes ribosome in train data)",
         fct_lvl = 1) %>%
  select(Recall, Precision, set, fct_lvl)

pr2 <- readr::read_csv(pr2f) %>%
  filter(label == "test") %>%
  mutate(set = "CFMS + APMS\n(ribosome PPIs removed from train data)",
         fct_lvl = 2) %>%
  select(Recall, Precision, set, fct_lvl)

pr3 <- readr::read_csv(pr3f) %>%
  filter(label == "test") %>%
  mutate(set = "CFMS + APMS\n(-ribosomes; +574123 missing PPI scores)",
         fct_lvl = 3) %>%
  select(Recall, Precision, set, fct_lvl) 
  
pr4 <- readr::read_csv(pr4f) %>%
  mutate(set = "CFMS + APMS\n(-ribosomes, +missing APMS, *NEW PIPELINE*)",
         fct_lvl = 4) %>%
  rename(Recall = recall, Precision = precision) %>% 
  select(Recall, Precision, set, fct_lvl)

all_pr <- rbind(pr1, pr2, pr3, pr4)

p <- ggplot(all_pr, aes(x = Recall, y = Precision, color = set, group = fct_rev(fct_reorder(set, fct_lvl)))) +
  geom_line(size = 1.5) +
  scale_color_manual(values = pal, name = "") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
p

ggsave("ppi_ml/figures/model_performance_w_fixedAPMS.png", p, device = "png", 
       width = 8, height = 10.5, units = "in")
