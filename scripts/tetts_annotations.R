library(tidyverse)
library(janitor)

`%!in%` = Negate(`%in%`)

# format annotations

leca_annots <- read_csv('leca/ppi_ml/annotations/leca_eunog_annots_complete.030721.csv') %>%
  select(-matches("*old*"), -matches("go_*"), -matches("*cmplx*"))

nog_mapping <- read_tsv('leca/ppi_ml/annotations/all.euNOG.diamond.mapping.2759') %>%
  filter(str_detect(ProteinID, 'TETTS')) %>% 
  mutate(id_fmt = str_extract(ProteinID,'(?<=\\|)(.*)(?=\\|)'))

tetts_uniprot <- read_tsv('/project/rmcox/misc/ty.tetts/annotations/tetts_uniprot.tsv') %>%
  clean_names() 

tetts_nog_annots <- tetts_uniprot %>%
  left_join(nog_mapping, by=c("entry"="id_fmt")) %>%
  select(ID, entry, everything()) %>%
  filter(!is.na(ID)) %>% 
  group_by(ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  left_join(leca_annots) 

tetts_annots <- tetts_uniprot %>%
  rename(ID = entry) %>% 
  filter(ID %!in% nog_mapping$id_fmt) %>%
  mutate(across(everything(), as.character)) %>% 
  bind_rows(tetts_nog_annots)

write_csv(tetts_annots, '/project/rmcox/misc/ty.tetts/annotations/tetts_complete_annots.csv')


# join annotations to mass spec data

tetts_ms <- read_tsv('/project/rmcox/misc/ty.tetts/results/clustered/tetts_cfms_normalized.clst.euclidean.average.elut') %>%
  select(-sparkline)

annot_cols <- names(tetts_annots)[! names(tetts_annots) %in% c('ID')]
data_cols <- names(tetts_ms)[! names(tetts_ms) %in% c('orthogroup')]
ordered_cols <- c("orthogroup", annot_cols, data_cols)

tetts_ms_annotated <- tetts_ms %>%
  left_join(tetts_annots, by=c("orthogroup"="ID")) %>%
  select(all_of(ordered_cols), -entry_name)
  #select(orthogroup, sparkline, entry, matches("*name*"), matches("*_cc"), matches("*gene_ontology*"), everything(), -entry_name, -ProteinID)

write_csv(tetts_ms_annotated, '/project/rmcox/misc/ty.tetts/results/clustered/tetts_cfms_normalized.clst.euclidean.average_annotated.csv')

  
