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
pr5f <- "ppi_ml/results/ppi_predict/all_feats/"


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

p <- ggplot(all_pr, aes(x = Recall, y = Precision, color = set, group = fct_reorder(set, fct_lvl))) +
  geom_line(size = 1.5) +
  scale_color_manual(values = pal, name = "") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow=2, byrow=TRUE, reverse=TRUE)) +
  geom_hline(yintercept = 0.85, linetype='dotted') +
  annotate("text", x = 0.85, y = 0.85, label = "15% FDR", vjust = -0.5)
p

ggsave("ppi_ml/figures/model_performance_new_pipeline.png", p, device = "png", 
       width = 12, height = 8, units = "in")


pr_5sel_f <- "ppi_ml/results/ppi_predict/5_feats/testset_precision_recall.csv"
pr_10sel_f <- "ppi_ml/results/ppi_predict/10_feats/testset_precision_recall.csv"
pr_20sel_f <- "ppi_ml/results/ppi_predict/20_feats/testset_precision_recall.csv"
pr_90sel_f <- "ppi_ml/results/ppi_predict/testset_precision_recall.csv"
pr_full_f <- "ppi_ml/results/ppi_predict/all_feats/testset_precision_recall.csv"


pr_5sel <- readr::read_csv(pr_5sel_f) %>%
  mutate(set = "5 features")
pr_10sel <- readr::read_csv(pr_10sel_f) %>%
  mutate(set = "10 features")
pr_20sel <- readr::read_csv(pr_20sel_f) %>%
  mutate(set = "20 features")
pr_91sel <- readr::read_csv(pr_90sel_f) %>%
  mutate(set = "91 features")
pr_full <- readr::read_csv(pr_full_f) %>%
  mutate(set = "663 features")

pr <- rbind(pr_5sel, pr_10sel, pr_20sel, pr_91sel, pr_full)
p <- ggplot(pr, aes(x = recall, y = precision, color = set)) +
  geom_line(size = 1, alpha = 0.75) +
  scale_color_manual(values = pal, name = "") +
  theme(legend.position = "bottom") +
  #guides(color = guide_legend(nrow=2, byrow=TRUE, reverse=TRUE)) +
  geom_hline(yintercept = 0.90, linetype='dotted') +
  annotate("text", x = 0.90, y = 0.90, label = "10% FDR", vjust = -0.5) #+
  #facet_wrap(~set, ncol=1)
p

compare <- pr %>%
  select(ID, ppi_score, set) %>%
  pivot_wider(names_from = set, values_from = ppi_score) %>%
  janitor::clean_names()

ggplot(compare, aes(x = x91_features, y = x10_features)) +
  geom_point(alpha = 0.5)
