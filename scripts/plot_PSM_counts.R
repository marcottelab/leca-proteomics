library(tidyverse)
library(network)
library(igraph)
library(tidygraph)
library(ggraph)
theme_set(theme_bw(base_size = 12))

pal_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
             "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
             "#7E6148", "#B09C85")
setwd("/stor/work/Marcotte/project/rmcox")

# ----------------------
# args
# ----------------------

ogs = c("KOG2386")
subset = c("human", "strpu", "xenla", "mouse", "chlre", "tetts", "brart", "euggr",
           "tryb2", "phatc")
pep_file = "leca/ppi_ml/data/cfms/pep_assign_posthoc_summarized.csv"
outdir = "leca/ppi_ml/figures/posthoc_peps"

# ----------------------
# global stuff
# ----------------------

# species names & clades
species <- read_csv('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/speciesinfo_clades.csv') %>%
  filter(tax_group == "eukaryota") %>% 
  select(code, species_name, clade) %>%
  mutate(code = tolower(code))
# ordered codes based on phylogeny
species_order <- read_delim('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/euk_codes_ordered_phylo.txt',
                            delim='\n', col_names=FALSE) %>% pull(X1)
# clade order
clade_order <- c('Amorphea','Excavate','TSAR','Archaeplastida')

library(tidyverse)
theme_set(cowplot::theme_cowplot())
pal_clades <- c("#E64B35", "#4DBBD5", "#3C5488", "#00A087")
rg <- c("#DC0000", "#979797")

pub_pal <- c("#E64B35", "#3578E6", "#4BD2AC", "#665899", "#E6A435","#8A8A95", "#353540")


# ----------------------
# script
# ----------------------

df <- read_csv(pep_file) %>%
  filter(orthogroup %in% ogs, !grepl("Fern", uniprot_id), 
         #species %in% subset,
         status == "unique", count > 50) %>%
  group_by(orthogroup) %>%
  arrange(desc(orthogroup), desc(count)) %>%
  left_join(species, by=c("species"="code"))

# lock in other var orders
df$clade <- factor(df$clade, levels = clade_order)
df$species <- factor(df$species, levels = species_order)

df %>%
  #filter(clade %in% c("Excavate", "TSAR")) %>%
  filter(clade == "Amorphea") %>% 
  ggplot(aes(x = gene_name, y = count, fill = clade)) +
    geom_col(width=0.75) +
    scale_fill_manual(values = pal_clades) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=8.5),
          strip.background = element_rect(color="#454545", 
                                          fill="#454545", 
                                          size=1.5, linetype="solid"),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "black", linetype = "dotted"),
          strip.text.x = element_text(color = "white", size = 8),
          strip.text.y = element_text(color = "white"),
          legend.title=element_blank(),
          legend.position = "top",
          legend.margin=margin(t=5),
          plot.margin = unit(c(0,0.5,0.5,0), "cm")) +
    labs(x = ogs, y = "Peptide spectral matches (unique)", fill = "") +
    #facet_wrap(~species, scales = "free_x", ncol=4)
    #facet_grid(species ~ orthogroup, scales = "free_x")
    facet_grid(~species, scales = "free_x", space = "free", switch = "x")

