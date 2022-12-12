infile = 'leca/ppi_ml/results/elutions/counts/cutoff_counts.csv'
chosen_cutoff = 150


cutoff_df = read_csv(infile)

max_ids = max(cutoff_df$num_ids)
max_psms = max(cutoff_df$cutoff)

final_ids <- cutoff_df %>%
  filter(cutoff == chosen_cutoff) %>%
  pull()

custom_label = paste(final_ids, 'proteins with >=', 
                     chosen_cutoff, 'PSMs',
                     sep = ' ')

theme_set(theme_bw(base_size=12))
p <- ggplot(cutoff_df, aes(x = cutoff, y = num_ids)) +
  geom_line() +
  geom_vline(xintercept = cutoff, linetype = "dashed") +
  scale_y_continuous(limits = c(0,max_ids), 
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(0,max_psms), 
                     expand = c(0, 0)) +
  annotate("text", x=chosen_cutoff+25, 
           y=(max_ids*.55), label=custom_label, 
           angle=270) +
  labs(x = "# PSMs", y = "# Proteins")
p

p %>% ggsave("leca/ppi_ml/figures/psm_thresholds.png", .,
         device = "png", width = 4, height = 4, 
         units = "in")
