library(tidyverse)

auc_file = 'leca/human_disease/network_propagation/results/leca_diseases_leave1out_auroc_top20_w-fu_sorted.csv'
annot_file = 'leca/ppi_ml/annotations/leca_eunog_annots_complete.030721.csv'

auc <- read_csv(auc_file)
annots <- read_csv(annot_file)

auc_fmt <- auc %>%
  select(disease, top_hits) %>%
  mutate(top_hits = strsplit(as.character(top_hits), ";")) %>% 
  unnest(top_hits) %>%
  left_join(annots, by=c("top_hits"="ID")) %>%
  select(-matches('old*'), -matches('*cmplx*'), -matches('*_up_ids*'))

auc_fmt$human_entry_counts <- sapply(strsplit(auc_fmt$human_entry,','), data.table::uniqueN)

write_csv(auc_fmt, 'leca/human_disease/network_propagation/results/top20_hits_annotated_w-fu.csv')

theme_set(theme_bw(base_size = 12))
pal <- c("#56B4E9", "#F0E442", "#009E24", "#E69F00", 
         "#FF0000", "#979797", "#5530AA", "#1E1E1E")
pal_lines <- c("#0072B2", "#E69F00")
auc %>%
  pivot_longer(cols = leave1out_auc:random_auc, names_to = "auc") %>%
  mutate(auc_label = ifelse(auc == "random_auc", "Random", "Back propagation")) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, alpha = 0.75, color = "black",
                 aes(fill = auc_label)) +
  scale_fill_manual(values = pal) +
  #geom_density(aes(color = auc_label)) +
  #geom_rug(aes(x=value, y=0, color=auc_label), position = position_jitter(height = 0)) +
  scale_color_manual(values = pal_lines) +
  guides(color = guide_legend(reverse=TRUE)) +
  labs(x = "area under ROC", y = "number of diseases (total = 109)") +
  theme(legend.title = element_blank(),
        legend.position = "top") -> p
p

p %>% ggsave("leca/human_disease/figures/network_prop_auc_w-fu_sorted.png", ., device = "png", width = 5, height = 3, units = "in")
