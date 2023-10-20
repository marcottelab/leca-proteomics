library(tidyverse)
library(janitor)
library(patchwork)
theme_set(theme_bw(base_size = 12))

palette_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85", "#42445A", "#2A303D",
                 "#3F404F", "#C03F2C")

workdir <- "/stor/work/Marcotte/project/rmcox/leca/"

plot_split_labels <- read_csv("dollo_qfo/annotations/eggnog/cog_slc_labels.csv") %>%
  select(category, plot_division)
leca_annots <- read_csv("dollo_qfo/annotations/leca_e5_annots.csv") %>%
  mutate(category_description = ifelse(category_description == 'Post-translational modification/protein turnover/chaperone functions',
                                       'Post-translational modification,\nprotein folding/turnover',
                                       category_description),
         category_description = str_replace(category_description, " and ", "/"),
         category_description = str_replace(category_description, "trafficing", "trafficking"),
         category_description = str_replace(category_description, "metabolis", "metabolism"),
         category_description = str_replace(category_description, "metabolismm", "metabolism"),
         category_description = str_replace(category_description, "envelop", "envelope"),
         category_description = str_replace(category_description, "Chromatin Structure/dynamics", "DNA replication/repair,\nchromatin dynamics"), 
         category_description = str_replace(category_description, "Replication/repair", "DNA replication/repair,\nchromatin dynamics"),
         category_description = str_replace(category_description, "Cytoskeleton", "Cytoskeleton/cell motility"),
         category_description = str_replace(category_description, "Cell motility", "Cytoskeleton/cell motility"),
         category_description = tolower(category_description),
         category_description = paste(toupper(substring(category_description, 1, 1)), substring(category_description, 2), sep = ""),
         category_description = str_replace(category_description, "Rna", "RNA"),
         category_description = str_replace(category_description, "Tranlsation", "Translation"), 
         category_description = str_replace(category_description, "Nucleotide metabolism/transport", "Nucleotide and coenzyme metabolism/transport"),
         category_description = str_replace(category_description, "Coenzyme metabolism", "Nucleotide and coenzyme metabolism/transport"),
         category_description = str_replace(category_description, "Secondary structure", "Secondary metabolite synthesis/transport")) %>%
  filter(category != "Y") %>% 
  left_join(plot_split_labels)

leca_annots %>% 
  ggplot(aes(x = fct_infreq(category_description), fill = human_status)) +
  geom_bar() +
  scale_fill_manual(values = c("#F39B7F","#42445A")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.7),
        legend.background = element_rect(colour = "#2A303D",linetype='solid')) +
  labs(y = '# of LECA OGs') -> p_all
p_all

p_all %>% ggsave("dollo_qfo/figures/leca_categories_101723.png", ., device = "png", width = 9.5, height = 4, units = "in")
p_all %>% ggsave("dollo_qfo/figures/leca_categories_101723.pdf", ., device = "pdf", width = 9.5, height = 4, units = "in")

# ------ top ------ 
leca_annots %>%
  filter(category_description != "Function unknown", plot_division=="T") %>% 
  ggplot(aes(x = fct_rev(fct_infreq(category_description)), fill = human_status)) +
  geom_bar() +
  scale_fill_manual(values = c("#F39B7F","#42445A")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                                   face="bold", size=10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "#2A303D",linetype='solid'),
        strip.background = element_blank(), strip.text = element_blank()) +
  labs(y = '# LECA OGs') +
  facet_grid(~plot_division, scales="free", space="free") -> p_top
p_top

p_top %>% ggsave("dollo_qfo/figures/leca_categories_top_101723.png", ., device = "png", width = 10, height = 5, units = "in")
p_top %>% ggsave("dollo_qfo/figures/leca_categories_top_101723.pdf", ., device = "pdf", width = 10, height = 5, units = "in")

# ------ bottom ------ 
leca_annots %>%
  filter(category_description != "Function unknown", 
         plot_division %in% c("N", "B")) %>% 
  ggplot(aes(x = fct_infreq(category_description), fill = human_status)) +
  geom_bar() +
  scale_fill_manual(values = c("#F39B7F","#42445A")) +
  theme(axis.text.x = element_text(angle = 45, hjust=0, vjust=0, 
                                   face="bold", size=10),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "#2A303D",linetype='solid'),
        strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_discrete(position = "top") +
  #facet_grid(~plot_division, scales="free", space="free") +
  labs(y = '# LECA OGs') -> p_bot
p_bot

p_bot %>% ggsave("dollo_qfo/figures/leca_categories_bot_101723.png", ., device = "png", width = 10, height = 5, units = "in")
p_bot %>% ggsave("dollo_qfo/figures/leca_categories_bot_101723.pdf", ., device = "pdf", width = 10, height = 5, units = "in")

# ------ unknown ------ 
leca_annots %>%
  filter(category_description == "Function unknown") %>% 
  ggplot(aes(x = fct_rev(fct_infreq(category_description)), fill = human_status)) +
  geom_bar() +
  scale_fill_manual(values = c("#F39B7F","#42445A")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "#2A303D",linetype='solid'),
        strip.background = element_blank(), strip.text = element_blank()) +
  labs(y = '# LECA OGs') +
  facet_grid(~plot_division, scales="free", space="free") -> p_unk
p_unk

p_unk %>% ggsave("dollo_qfo/figures/leca_categories_unk_101723.png", ., device = "png", width = 1.5, height = 5, units = "in")
p_unk %>% ggsave("dollo_qfo/figures/leca_categories_unk_101723.pdf", ., device = "pdf", width = 1.5, height = 5, units = "in")
