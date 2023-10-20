library(tidyverse)
theme_set(cowplot::theme_cowplot())
pal_clades <- c("#E64B35", "#4DBBD5", "#3C5488", "#00A087")

draft_labels = read_csv("leca/ppi_ml/results/walktrap/LinearSVC_100feats_fdr10_4steps_nochloro_explode_algo_labels_091923_fixed.csv")

cuts2eval = c("cut_796_algo_label", "cut_1202_algo_label","cut_1608_algo_label","cut_2014_algo_label")
df <- draft_labels %>%
  pivot_longer(cut_557_algo_label:cut_2014_algo_label, names_to='cut', values_to='label') %>% 
  group_by(cut, label) %>%
  filter(cut %in% cuts2eval) %>% 
  mutate(categ = case_when(!label %in% c("Unclustered", "Large heterogeneous complex", "Uncharacterized complex") ~ "Known complex label",
                           TRUE ~ label))
df

df$cut <- factor(df$cut, levels = cuts2eval)
ggplot(df, aes(x = cut, fill = categ)) +
  geom_bar() +
  scale_fill_manual(values = pal_clades)
