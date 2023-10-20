library(tidyverse)
library(janitor)
theme_set(theme_bw(base_size = 12))

palette_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85", "#42445A", "#2A303D",
                 "#3F404F", "#C03F2C")

workdir <- "/stor/work/Marcotte/project/rmcox/leca/"
data_path <- paste0(workdir,'localization_ml/results/leca_prop/')
files <- dir(data_path, pattern = "*.csv")


# ------------------------------------------------
# meta data
# ------------------------------------------------
# total # of LECA OGs: 10092
# total # of human LECA GO term assignments: 29674
# total # of arath LECA GO term assignments: 16133
# total # of yeast LECA GO term assignments: 7447
# total # of tetts LECA GO term assignments: 1493
# total # of tryb2 LECA GO term assignments: 7203
# ------------------------------------------------

# read in data
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>% 
  reduce(rbind) %>%
  clean_names()

# combine some labels
cilia = c("cilium", "motile cilium")
combined_cilium <- data %>%
  filter(go_name %in% cilia) %>% 
  group_by(leca, species) %>%
  summarize(count = sum(count)) %>%
  mutate(go_name = "cilium")

transport_vesicles = c("COPII-coated ER to Golgi transport vesicle", 
                       "COPI-coated vesicle", "clathrin-coated vesicle")
combined_tv <- data %>%
  filter(go_name %in% transport_vesicles) %>% 
  group_by(leca, species) %>%
  summarize(count = sum(count)) %>%
  mutate(go_name = "transport vesicle")

digestive_vesicles = c("lysosome", "phagocytic vesicle", "autophagosome")
combined_dv <- data %>%
  filter(go_name %in% digestive_vesicles) %>% 
  group_by(leca, species) %>%
  summarize(count = sum(count)) %>%
  mutate(go_name = "digestive vesicle")

# format data
species_labels = c("S. cerevisiae (n=7447)", "H. sapiens (n=29674)",
                   "T. brucei (n=7203)", "T. thermophila (n=1493)", 
                   "A. thaliana (n=16133)")

to_remove = c("melanosome", "phagocytic cup", "contractile vacuole", 
              "chromaffin granule", "cell tip", cilia, 
              transport_vesicles, digestive_vesicles)

plot_data <- data %>%
  filter(!go_name %in% to_remove) %>%
  bind_rows(combined_cilium, combined_tv, combined_dv) %>%
  mutate(label = ifelse(leca==TRUE, "LECA or older", "Younger than LECA"),
         species_fmt = case_when(species=="yeast" ~ species_labels[1],
                                 species=="human" ~ species_labels[2],
                                 species=="tryb2" ~ species_labels[3],
                                 species=="tetts" ~ species_labels[4],
                                 species=="arath" ~ species_labels[5]),
         go_name = ifelse(go_name=="perinuclear region of cytoplasm", "perinuclear region",
                          go_name),
         go_name = paste(toupper(substring(go_name, 1, 1)), substring(go_name, 2), sep = "")
  ) %>%
  select(go_name, species_fmt, label, count) %>% 
  group_by(go_name, species_fmt) %>%
  mutate(perc = (round(count / sum(count), 3))*100,
         perc_label = paste0(perc, '%')) %>%
  ungroup() %>%
  mutate(summary = ifelse(label=="LECA or older", paste0(perc_label, " (", count, "/", (round(count/(perc/100), 0)), ")"), NA))

# refactor variables for plot ordering
plot_data$species_fmt = factor(plot_data$species_fmt, levels = species_labels)
desc_slc <- plot_data %>%
  filter(species_fmt=="H. sapiens (n=29674)", label=="LECA or older") %>%
  arrange(perc) %>%
  pull(go_name)
plot_data$go_name = factor(plot_data$go_name, levels = desc_slc)

# plot results
plot_data %>%
  ggplot(aes(x = go_name, y = count)) +
  geom_col(aes(fill = factor(label, levels=c("Younger than LECA","LECA or older"))),
               position = "fill") +
  #geom_label(aes(x = go_name, group = label, label = summary), position=position_fill(reverse=TRUE, vjust=0.2), size = 1.5) +
  scale_fill_manual(values = c("#91D1C2","#42445A"),
                    name = "Estimated gene age") +
  scale_y_continuous(breaks=seq(0,1,0.5)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~species_fmt, nrow = 1) +
  coord_flip() +
  theme(legend.position = "top",
        legend.key.width = unit(0.5, "in"),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold"),
        strip.text = element_text(size=9.25), # facet text size
        axis.title.y = element_blank()
        ) +
  ylab("Proportion (*genes can be assigned to multiple GO terms)") -> p
p

p %>% ggsave("localization_ml/figures/leca_localization.png", ., device = "png", width = 11.4, height = 4.5, units = "in")
p %>% ggsave("localization_ml/figures/leca_localization.pdf", ., device = "pdf", width = 11.4, height = 4.5, units = "in")
  
 
