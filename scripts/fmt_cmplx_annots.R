library(tidyverse)

leca_nog_map <- read_tsv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/all.euNOG.diamond.mapping.2759") %>%
  mutate(uniprot_ID = str_extract(ProteinID,'(?<=\\|)(.*)(?=\\|)')) %>%
  select(-ProteinID) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', ')))

cp_annots <- read_csv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_comportal.csv") %>%
  janitor::clean_names() %>%
  select(ids_fmt, number_complex_ac, recommended_name, 
         aliases_for_complex, complex_assembly) %>%
  mutate(ids = strsplit(as.character(ids_fmt), " ")) %>% 
  unnest(ids) %>%
  unique() %>%
  rename(uniprot_ID = ids, 
         cp_cmplx_name = recommended_name, 
         cp_cmplx_synonym = aliases_for_complex, 
         cp_cmplx_members = ids_fmt,
         cp_cmplx_assembly = complex_assembly) %>%
  select(uniprot_ID, cp_cmplx_name, cp_cmplx_synonym, 
         cp_cmplx_members, cp_cmplx_assembly) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  left_join(leca_nog_map) %>%
  filter(!is.na(ID)) %>%
  mutate(ID = strsplit(as.character(ID), ", ")) %>% 
  unnest(ID) %>%
  group_by(ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  select(ID, everything()) %>%
  rename(cp_up_ids = uniprot_ID)
  
  
corum_annots <- read_tsv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_corum.txt") %>%
  janitor::clean_names() %>%
  rename(corum_cmplx_name = complex_name, 
         corum_cmplx_members = subunits_uni_prot_i_ds, 
         corum_cmplx_synonym = synonyms, 
         corum_cmplx_comment = complex_comment,
         corum_cmplx_disease = disease_comment) %>%
  mutate(uniprot_ID = strsplit(as.character(corum_cmplx_members), ";")) %>%
  unnest(uniprot_ID) %>% 
  unique() %>% 
  select(uniprot_ID, matches('corum*')) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  left_join(leca_nog_map) %>%
  filter(!is.na(ID)) %>%
  mutate(ID = strsplit(as.character(ID), ", ")) %>% 
  unnest(ID) %>%
  group_by(ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  select(ID, everything()) %>%
  rename(corum_up_ids = uniprot_ID)

go_annots <- read_tsv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_GO_human.tsv") %>% 
  janitor::clean_names() %>%
  rename(uniprot_ID = gene_product_id, go_cmplx_name = go_name,
         go_gene_name = symbol) %>%
  select(uniprot_ID, go_cmplx_name, go_gene_name) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  left_join(leca_nog_map) %>%
  filter(!is.na(ID)) %>%
  mutate(ID = strsplit(as.character(ID), ", ")) %>% 
  unnest(ID) %>%
  group_by(ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  select(ID, everything()) %>%
  rename(go_up_ids = uniprot_ID)

all_cmplx_annots <- go_annots %>%
  full_join(corum_annots, by=c("ID")) %>%
  full_join(cp_annots, by=c("ID")) %>%
  select(ID, go_gene_name, matches('*cmplx_name*'), matches('*cmplx_synonym*'),
         matches('*up_ids*'), everything())

write_csv(all_cmplx_annots, "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/all_cmplx_annots.csv")




### JOIN COMPLEX ANNOTS ON UNIPROT ID

go_annots <- read_tsv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_GO_human.tsv") %>% 
  janitor::clean_names() %>%
  rename(uniprot_ID = gene_product_id, go_cmplx_name = go_name,
         go_gene_name = symbol) %>%
  select(uniprot_ID, go_cmplx_name, go_gene_name) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', ')))

corum_annots <- read_tsv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_corum.txt") %>%
  janitor::clean_names() %>%
  rename(corum_cmplx_name = complex_name, 
         corum_cmplx_members = subunits_uni_prot_i_ds, 
         corum_cmplx_synonym = synonyms, 
         corum_cmplx_comment = complex_comment,
         corum_cmplx_disease = disease_comment) %>%
  mutate(uniprot_ID = strsplit(as.character(corum_cmplx_members), ";")) %>%
  unnest(uniprot_ID) %>% 
  unique() %>% 
  select(uniprot_ID, matches('corum*')) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', ')))

cp_annots <- read_csv("/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_comportal.csv") %>%
  janitor::clean_names() %>%
  select(ids_fmt, number_complex_ac, recommended_name, 
         aliases_for_complex, complex_assembly) %>%
  mutate(ids = strsplit(as.character(ids_fmt), " ")) %>% 
  unnest(ids) %>%
  unique() %>%
  rename(uniprot_ID = ids, 
         cp_cmplx_name = recommended_name, 
         cp_cmplx_synonym = aliases_for_complex, 
         cp_cmplx_members = ids_fmt,
         cp_cmplx_assembly = complex_assembly) %>%
  select(uniprot_ID, cp_cmplx_name, cp_cmplx_synonym, 
         cp_cmplx_members, cp_cmplx_assembly) %>%
  group_by(uniprot_ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', ')))

all_cmplx_annots <- go_annots %>%
  full_join(corum_annots, by=c("uniprot_ID")) %>%
  full_join(cp_annots, by=c("uniprot_ID"))

write_csv(all_cmplx_annots, "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/cmplx/cmplx_annots_all_upIDs.csv")

