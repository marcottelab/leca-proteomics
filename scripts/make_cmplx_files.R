setwd("/project/rmcox/LECA/ms/gold_stds/")
`%!in%` = Negate(`%in%`)

# how to map TAIR IDs to eggNOG groups?
# both of these mapping files are highly redundant
arath_idmapping <- read_tsv("tair_up_codes.txt") %>%
 select(uniprot_id, name) %>%
 mutate(name = str_replace(name, "\\.\\d", "")) %>%
  unique()

arath_idmapping %>% 
  group_by(name) %>% 
  tally() %>%
  arrange(desc(n)) %>%
  filter(name %in% tair$locus_name)

tair_idmapping <- read_tsv("TAIR2UniprotMapping.txt", col_names = FALSE) 

tair_idmapping %>%
  group_by(X1) %>% 
  tally() %>%
  arrange(desc(n))

# list of proteomes already mappedd to eggNOG groups
species <- read_csv("gold_std_info.csv")

# gold std protein complexes for many species from Complex Portal 
cp_dir <- "/project/rmcox/LECA/ms/gold_stds/comportal"
cp_files <- dir(cp_dir, pattern = "*[[:digit:]].gold.cmplx.txt")
comportal <- cp_files %>%
  map(~ read_tsv(file.path(cp_dir, .), col_names = FALSE)) %>% 
  reduce(rbind)

  comportal %>%
    write_delim("comportal/all.gold.cmplx.txt")

# gold std mammalian complexes from CORUM
corum_all <- read_tsv("corum/allComplexes.txt") %>%
  janitor::clean_names()

  corum_all %>%
    select(subunits_uni_prot_i_ds) %>%
    mutate(ids_fmt = str_replace_all(subunits_uni_prot_i_ds, ";", " ")) %>% 
    select(ids_fmt) %>%
    write_delim("corum/corum.cmplx.txt", col_names = FALSE)

# gold std Arabidopsis pairwaise interactions from TAIR
tair <- read_tsv("tair/TairProteinInteraction.20090527.txt") %>%
  janitor::clean_names() 

  tair_types <- tair %>% pull(evidence_description) %>% unique()
  tair_types_bad <- c(tair_types[5], tair_types[12], tair_types[21], tair_types[25])
  
  tair_fmt <- tair %>%
    filter(evidence_description %!in% tair_types_bad) %>%
    select(locus_name, interactor_locus_name) %>%
    left_join(arath_idmapping, by = c("locus_name" = "name")) %>%
    rename(p1 = uniprot_id) %>% 
    left_join(arath_idmapping, by = c("interactor_locus_name" = "name")) %>%
    rename(p2 = uniprot_id) %>% 
    unique() %>%
    mutate(ids_fmt = paste(p1, p2, " "))
  
  tair_fmt %>%
    select(ids_fmt) %>%
    write_delim("tair/tair.cmplx.txt")

# manually curated Arabidopsis protein complexes by Claire
claire <- read_tsv("claire/arath_manual.txt")

  claire %>%
    select(complex) %>%
    write_delim("claire/claire.cmplx.txt")
  

# combine all protein complexes
comportal_nog <- read_delim("comportal/all.euNOG.gold.cmplx.txt", delim = "\n", col_names = F)
corum_nog <- read_delim("corum/corum.euNOG.gold.cmplx.txt", delim = "\n", col_names = F)
claire_nog <- read_delim("claire/claire.euNOG.gold.cmplx.txt", delim = "\n", skip = 1, col_names = F)
tair_nog <- read_delim("tair/tair.euNOG.gold.cmplx.txt",  delim = "\n", col_names = F)

all_cmplx <- rbind(comportal_nog, corum_nog, claire_nog, tair_nog)
write_delim(all_cmplx, "all.gold.cmplx.txt", col_names = F)

# filter out ribosomal complexes
ribokogs <- read_delim("ribosomal_kogs.txt", delim = "\n", col_names = F)

not_ribo <- all_cmplx %>%
  mutate(ribosomal = case_when(grepl(paste(ribokogs$X1, collapse = "|"), X1) ~ "TRUE",
                               TRUE ~ "FALSE")) %>%
  filter(ribosomal == "FALSE") %>%
  select(X1) %>%
  write_delim("all.gold.cmplx.noRibos.txt", col_names = F, delim = "\n")
  


