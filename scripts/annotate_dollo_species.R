library(ggtree)
library(ape)
library(treeio)
library(tidytree)
library(janitor)
library(stringi)
#library(phlyobase)

# recreate meta data
# should be 156 total species
count_file = '../leca/dollo_qfo/results/qfo_count_tables.tsv'
tax2up_file = '../leca/dollo_qfo/data/proteomes/Added_Species/spec2taxid_mapping.txt'
dollo_data <- read_tsv(count_file)
tax2up <- read_delim(tax2up_file, delim=" ", col_names=F) %>%
  rename(species_code = X1, tax_id = X2)

dollo_species <- dollo_data %>%
  pivot_longer(!family, names_to = "species_code", values_to = "presence") %>%
  select(species_code) %>%
  unique()

all_species <- dollo_species %>%
  left_join(tax2up)

# join recovered meta data
clean_meta <- metadata %>%
  clean_names() %>%
  mutate(tax_group = str_extract(taxonomic_lineage, "[^, ]+")) %>%
  select(-taxonomic_lineage) %>%
  unique() %>%
  mutate(organism = str_replace(organism, " / ", ""),
         species_name = str_extract(organism, "^[^(]*")) %>%
  filter(taxon_mnemonic != "null")

write_csv(clean_meta, '../leca/dollo_qfo/annotations/species_info_100223.csv')

clean_meta %>% select(species_name) %>% write_delim('../leca/dollo_qfo/data/full_species_list.txt', col_names=F)

# tree info
nhx_file = '../leca/dollo_qfo/data/species_tree/speciestree.nhx'
nwk_file = '../leca/dollo_qfo/data/species_tree/speciestree.nwk'
annot_file = '../leca/dollo_qfo/annotations/species_info_100223.csv'

tree <- read.tree(nwk_file)
tree_df <- as_tibble(tree)

annot_df <- read_csv(annot_file)

tree_annots <- tree_df %>%
  mutate(tax_id = as.integer(str_extract(label, "(?<=__)\\d+$"))) %>%
  left_join(annot_df, by=c("tax_id"="organism_id")) %>%
  filter(!is.na(proteome_id))

all_species %>%
  filter(!species_code %in% tree_annots$taxon_mnemonic)
