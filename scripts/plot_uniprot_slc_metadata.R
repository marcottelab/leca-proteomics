library(tidyverse)
library(janitor)
library(patchwork)
library(cowplot)
library(scales)
theme_set(theme_bw(base_size = 12))

palette_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85", "#42445A", "#2A303D",
                 "#3F404F", "#C03F2C")

workdir <- "/stor/work/Marcotte/project/rmcox/leca/"

## -------------------------------------------------------------
## functions
## -------------------------------------------------------------

fmt_slc_data <- function(file, colors) {
  
  slc_meta <- str_extract(file, "(?<=/annotations/).*") %>%
    str_replace(".tsv", "") %>% 
    str_split("_SL-", n=2, simplify=T)
  
  slc = str_replace_all(slc_meta[1], "_", " ")
  slc_id = paste0("SL-", slc_meta[2])
  
  df <- read_tsv(file) %>%
    clean_names() %>%
    filter(!is.na(entry_name)) %>% 
    select(entry, entry_name, organism_id, protein_families, gene_ontology_go,
           subcellular_location_cc, taxonomic_lineage) %>% 
    mutate(taxonomic_lineage = str_replace_all(taxonomic_lineage, 
                                               " [\\(\\[].*?[\\)\\]]", "")) %>%
    separate(taxonomic_lineage, into=c("origin","branch","clade",
                                       "subclade_1","subclade_2"), sep=", ") %>%
    mutate(label = case_when(clade %in% amorphea ~ "Amorphea",
                             clade %in% archaeaplastida ~ "Archaeaplastida",
                             clade %in% tsar ~ "TSAR",
                             clade %in% excavata ~ "Excavata",
                             clade %in% cryptophyta ~ "Cryptophyta",
                             clade %in% haptophyta ~ "Haptophyta",
                             branch %in% viridae ~ "Viridae",
                             TRUE ~ branch)) %>%
    mutate(slc = slc, slc_id = slc_id) %>%
    left_join(colors)
  return(df)
}

## -------------------------------------------------------------
## labels
## -------------------------------------------------------------

# labels to remove
to_remove = c("Caveola", "Chromaffin granule", "Melanosome", "Porosome", 
              "Cell tip", "Cell cortex", "Contractile vacuole")

# labels to highlight
to_highlight = c("Phagosome", "Phagocytic cup", "Lipid droplet")

# labels to combine
cilium = c("Flagellum", "Cilium")
cargo_ves = c("COPI-coated vesicle", "COPII-coated vesicle", 
              "Clathrin-coated vesicle")
digestive_ves = c("Lysosome", "Phagosome", "Peroxisome")
podia = c("Pseudopodium", "Lamellipodium", "Filopodium")
cyto_granule = c("Cytolytic granule", "Stress granule",
                 "Cytoplasmic ribonucleoprotein granule")

## -------------------------------------------------------------
## data sets
## -------------------------------------------------------------

# read in SLC nog mapping data
mapping <- read_csv("localization_ml/results/uniprot_leca_slc.csv") %>%
  mutate(leca_label = ifelse(leca==1, "LECA or older", "Younger than LECA"))

# read in and format SLC annotation data
amorphea = c("Opisthokonta","Amoebozoa")
archaeaplastida = c("Viridiplantae","Rhodophyta","Glaucocystophyceae")
tsar = c("Sar")
excavata = c("Discoba","Metamonada")
cryptophyta = c("Cryptophyceae")
haptophyta = c("Haptista")
viridae = c("Duplodnaviria", "Naldaviricetes", "Riboviria", 
            "Varidnaviria", "unclassified bacterial viruses")

label = c("Amorphea", "Excavata", "TSAR", "Haptophyta", "Cryptophyta", "Archaeaplastida", 
          "Archaea", "Bacteria", "Viridae")
color = c("#E64B35", "#4DBBD5", "#3C5488", "#006471", "#91D1C2", "#00A087",
          "#E69F00", "#64657A", "black")
clade_colors = data.frame(label, color)

data_path = "localization_ml/data/uniprot_slc/annotations"
files = list.files(data_path, "*.tsv")
slc_df <- files %>%
  map(~ fmt_slc_data(file.path(data_path, .), clade_colors)) %>%
  reduce(rbind) %>%
  mutate(corrected_slc = ifelse(str_detect(tolower(subcellular_location_cc), 
                                           tolower(slc)), slc, NA)) # quality control slc assignments

# problem IDs...
error_count <- slc_df %>% 
  filter(is.na(corrected_slc)) %>%
  group_by(slc) %>%
  tally()

# calculate counts per SLC
slc_nogs <- slc_df %>% 
  left_join(mapping, by=c("entry", "slc")) 

slc_nog_counts <- slc_nogs %>%
  filter(!is.na(corrected_slc), !is.na(ID)) %>%
  select(-slc) %>% rename(slc = corrected_slc) %>% 
  mutate(slc = ifelse(slc %in% cargo_ves, "Transport vesicle", slc),
         slc = ifelse(slc %in% cilium, "Cilium", slc),
         slc = ifelse(slc %in% cyto_granule, "Cytoplasmic granule", slc),
         slc = ifelse(slc %in% podia, "Cell projection", slc)) %>% 
  filter(!slc %in% to_remove) %>%
  select(slc, ID) %>%
  unique() %>% 
  group_by(slc) %>% tally() %>%
  ungroup() %>%
  rename(n_eunog = n) %>% 
  arrange(desc(n_eunog))
  
# read in and format count data
data <- read_tsv("localization_ml/data/uniprot_slc/uniprot_slc_counts.tsv", 
                 col_names=F) %>%
  rename(count=X1, file=X2) %>%
  separate(file, into=c("slc","slc_id"), sep="_SL-", remove=F) %>%
  mutate(slc = str_replace_all(slc, "_", " "),
         slc_id = str_replace(slc_id, ".tsv", ""),
         slc_id = paste0("SL-", slc_id)) %>%
  left_join(error_count) %>%
  mutate(corrected_count = count - n) %>% 
  mutate(slc = ifelse(slc %in% cargo_ves, "Transport vesicle", slc),
         slc = ifelse(slc %in% cilium, "Cilium", slc),
         slc = ifelse(slc %in% cyto_granule, "Cytoplasmic granule", slc),
         slc = ifelse(slc %in% podia, "Cell projection", slc)) %>% 
  filter(!slc %in% to_remove) %>%
  group_by(slc) %>% summarize(total_count = sum(corrected_count)) %>%
  ungroup() %>%
  arrange(desc(total_count)) %>%
  left_join(slc_nog_counts)

# set category plot order
slc_levels <- pull(data, slc) %>% unique() %>% rev()
data$slc <- factor(data$slc, levels=slc_levels)

## -------------------------------------------------------------
## counts per SLC
## -------------------------------------------------------------

# lollipop plot
p1 <- ggplot(data, aes(x=slc)) +
  geom_segment(
    aes(x=slc, xend=slc, y=0, yend=total_count), 
    #color=ifelse(data$slc %in% to_highlight, "#91D1C2", "#42445A"),
    #size=ifelse(data$slc %in% to_highlight, 1.75, 1.75)
    color="grey", size=1.75
  ) +
  geom_point(aes(y=total_count, color="total"),
    #color=ifelse(data$slc %in% to_highlight, "#91D1C2", "#42445A"), 
    #size=ifelse(data$slc %in% to_highlight, 5, 5)
    size=5
  ) +
  geom_segment(
    aes(x=slc, xend=slc, y=0, yend=n_eunog), 
    color="#3F404F", size=1.75
  ) +
  geom_point(aes(y = n_eunog, color="nog"),
    #color=ifelse(data$slc %in% to_highlight, "#91D1C2", "#42445A"), 
    #size=ifelse(data$slc %in% to_highlight, 5, 5)
    size=5
  ) +
  scale_y_log10(labels = scales::comma) +
  theme_bw(base_size = 14) +
  coord_flip() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y = element_text(face="bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box="vertical", legend.margin=margin()
  ) +
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual(breaks=c("nog", "total"),
                     values=c(nog='#3F404F', total='grey'),
                     labels=c(nog='# unique euNOGs',
                              total='total # associated proteins\nin UniProt database (reviewed)'))
p1

# add annotation for highlighted SL terms
# data_highlight <- data %>%
#   filter(slc %in% to_highlight)
# hl_labels <- pull(data_highlight, slc)
# hl_counts <- pull(data_highlight, total_count)
# p1 <- p1 + annotate("text", x=hl_labels, y=hl_counts, 
#               label=paste("n = ", hl_counts, sep=""),
#               color="#91D1C2", size=4 , angle=0, fontface="bold", hjust=-0.3)
# p1

## -------------------------------------------------------------
## species per SLC
## -------------------------------------------------------------

# get counts per clade per slc
counts <- slc_df %>%
  select(-slc) %>% 
  rename(slc = corrected_slc) %>% 
  filter(!slc %in% to_remove) %>%
  #filter(!label %in% c("Bacteria", "Viridae")) %>% 
  group_by(label, slc, color) %>%
  mutate(slc = ifelse(slc %in% cargo_ves, "Transport vesicle", slc),
         slc = ifelse(slc %in% cilium, "Cilium", slc),
         slc = ifelse(slc %in% cyto_granule, "Cytoplasmic granule", slc),
         slc = ifelse(slc %in% podia, "Cell projection", slc)) %>%
  tally() %>%
  filter(!is.na(slc))

# get counts per clade per subphyla per slc
counts_granular <- slc_df %>%
  select(-slc) %>% 
  rename(slc = corrected_slc) %>% 
  filter(!slc %in% to_remove) %>%
  #filter(!label %in% c("Bacteria", "Viridae")) %>% 
  group_by(label, clade, slc, color) %>%
  mutate(slc = ifelse(slc %in% cargo_ves, "Transport vesicle", slc),
         slc = ifelse(slc %in% cilium, "Cilium", slc),
         slc = ifelse(slc %in% cyto_granule, "Cytoplasmic granule", slc),
         slc = ifelse(slc %in% podia, "Cell projection", slc)) %>%
  tally() %>%
  filter(!is.na(slc))

# get species counts per clade
species_counts <- slc_df %>%
  select(-slc) %>% 
  rename(slc = corrected_slc) %>% 
  filter(!slc %in% to_remove) %>%
  mutate(slc = ifelse(slc %in% cargo_ves, "Transport vesicle", slc),
         slc = ifelse(slc %in% cilium, "Cilium", slc),
         slc = ifelse(slc %in% cyto_granule, "Cytoplasmic granule", slc),
         slc = ifelse(slc %in% podia, "Cell projection", slc)) %>%
  select(organism_id, label, slc, color) %>%
  unique() %>%
  group_by(label, slc, color) %>%
  tally()

# set category plot order
slc_levels <- data %>% arrange(desc(total_count)) %>%
  pull(slc) %>% unique() %>% rev()
counts$slc <- factor(counts$slc, levels=slc_levels)
counts$label <- factor(counts$label, levels=label)
species_counts$slc <- factor(species_counts$slc, levels=slc_levels)
species_counts$label <- factor(species_counts$label, levels=label)

# dot plot
p2 <- counts %>%
  filter(label != "Viridae") %>% 
  ggplot(aes(x = label, y = slc, color = color)) +
  geom_point(shape = 1, color = "black", size=1.75, stroke=3) +
  geom_point(aes(color = color,
                 size=ifelse(n>0, 1.75, NA))) +
  theme_bw(base_size = 14) +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(t=-10)) +
  scale_size(guide = 'none') +
  labs(color = "")
p2

p2_spe <- species_counts %>%
  filter(label != "Viridae",
         !is.na(slc)) %>%
  mutate(n = ifelse(n>1000, 1000, n)) %>% 
  ggplot(aes(x = label, y = slc, color = color)) +
  #geom_point(shape = 1, color = "black", size=n, stroke=1) +
  geom_point(aes(color = color,
                 size=n)) +
  theme_bw(base_size = 14) +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 35, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(t=-10),
        ) +
  scale_size(range = c(1,5),
             limits = c(min(species_counts$n), 1000),
             breaks = c(1, 100, 1000),
             labels = c("1", "100", ">=1000"),
             name = "# unique species") +
  labs(color = "")
p2_spe

## -------------------------------------------------------------
## number LECA per SLC
## -------------------------------------------------------------
leca_e5_unknown <- read_csv("dollo_qfo/annotations/leca_e5_annots.csv") %>%
  filter(category_description=="Function Unknown")

mapping_fmt <- mapping %>% 
  mutate(slc = ifelse(slc %in% cargo_ves, "Transport vesicle", slc),
         slc = ifelse(slc %in% cilium, "Cilium", slc),
         slc = ifelse(slc %in% cyto_granule, "Cytoplasmic granule", slc),
         slc = ifelse(slc %in% podia, "Cell projection", slc)) %>% 
  filter(!slc %in% to_remove) %>%
  mutate(leca_label = ifelse(ID %in% leca_e5_unknown$euNOG, "LECA or older & 'function unknown'", leca_label)) %>%
  select(ID, leca_label, slc) %>%
  unique() %>%
  group_by(leca_label, slc) %>%
  tally() %>%
  arrange(desc(n))

mapping_fmt$slc <- factor(mapping_fmt$slc, levels=slc_levels)
mapping_fmt <- filter(mapping_fmt, !is.na(slc))

leca_counts <- mapping_fmt %>%
  filter(leca_label != "Younger than LECA") %>%
  group_by(slc) %>%
  summarize(total_leca = sum(n)) %>%
  mutate(leca_label = "LECA or older",
         leca_n = paste(total_leca, "LECA OGs", sep=" "))

# total euNOGs
mapping %>% 
  select(ID, leca_label) %>% 
  unique() %>% tally()

# total leca representation ~ 5312/10092
mapping %>% 
  select(ID, leca_label) %>% 
  filter(leca_label != "Younger than LECA") %>% 
  unique() %>% tally() 

# bar chart
p3 <- ggplot(mapping_fmt, aes(x=slc, y=n,
                              fill=fct_rev(leca_label))) +
  geom_col(width=0.7, position="fill", size=0.8, color="black") +
  # geom_text(data=leca_counts, aes(label=leca_n, x=slc),
  #           y=0.10,
  #           color = "white") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.box="vertical", legend.margin=margin()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("#C6C4C9","#F8AA81","#3F404F")) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 2)) +
  coord_flip() 
p3

## -------------------------------------------------------------
## combined panel
## -------------------------------------------------------------

p1_strip <- p1 + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank())
p2_strip <- p2_spe + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
p3_strip <- p3 + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())

panel <- (p1_strip + p2_strip + p3_strip) + 
  plot_layout(nrow=1, widths = c(1, 1, 1.2)) + 
  #plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "bottom")
panel

panel %>% ggsave("localization_ml/figures/leca_uniprot_counts_species_102223.png", ., device = "png", width = 12, height = 6, units = "in")
panel %>% ggsave("localization_ml/figures/leca_uniprot_counts_species_102223.pdf", ., device = "pdf", width = 12, height = 6, units = "in")

## -------------------------------------------------------------
## write out leca gene sets
## -------------------------------------------------------------

# slc_nogs <- slc_df %>% 
#   left_join(mapping, by=c("entry", "slc")) %>%
#   filter(!is.na(corrected_slc), !is.na(ID))
# leca <- filter(slc_nogs, leca==1)
# unique_slc <- leca %>% pull(corrected_slc) %>% unique()
# 
# for (i in unique_slc){
#   
#   print(paste0(i, ":"))
#   fname = str_replace_all(i, " ", "_")
#   outdir = "localization_ml/results/leca_slc_gene_sets/"
#   outname = paste0(outdir, fname, ".csv")
#   df <- leca %>%
#     filter(corrected_slc==i) %>%
#     write_csv(outname)
#   print(df)
#   
# }
