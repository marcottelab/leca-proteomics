library(tidyverse)
theme_set(theme_bw(base_size = 12))

# make AA <-> color mapping
colors = c("#4BD2AC","#E6A435", "#3578E6", "#FDE725", "#FA5FFA", "#665899", "#21918C")
aa = c("S,T,Q", "G", "I,V,L,F,A", "P", "R", "D,E", "H")

aa_df <- data.frame(colors, aa) %>%
  separate_rows(aa)

# read in and format seq data 
seqs <- read_delim('leca/ppi_ml/highlights/efhc2/efhc2_seqs.txt', delim=' ') %>%
  pivot_longer(`1`:`23`, names_to = "position", values_to = "aa") %>%
  left_join(aa_df)

# species order
spec_order <- c("HS", "MM", "GG", "XT", "DR", "CI", "CE", "DM", "TT", "CR", "TB")
seqs$sp <- factor(seqs$sp, levels = spec_order)

# generate plot
plot <- seqs %>% filter(!sp %in% c("GG","DR","CI")) %>% 
  ggplot(aes(x = as.numeric(position), y = fct_rev(sp),
             label = aa)) +
  geom_tile(aes(fill = colors), alpha = 0.9) +
  scale_fill_identity() +
  geom_text(fontface = "bold") +
  scale_x_discrete(breaks = c(1,23)) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )
plot

# save plot
prefix = "leca/human_disease/figures/efhc2_seqs"
h = 3.6; w = 5.2
plot %>% ggsave(paste0(prefix, ".png"), ., device = "png",
                width = w, height = h, units = "in")
plot %>% ggsave(paste0(prefix, ".pdf"), ., device = "pdf",
                width = w, height = h, units = "in")