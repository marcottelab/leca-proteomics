cp_annots <- read_csv("/stor/work/Marcotte/project/rmcox/LECA/ms/gold_stds/comportal/all.annotations.gold.cmplx.csv") %>%
  janitor::clean_names() %>%
  select(ids_fmt, number_complex_ac, recommended_name, 
         aliases_for_complex, complex_assembly) %>%
  mutate(ids = strsplit(as.character(ids_fmt), " ")) %>% 
  unnest(ids) %>%
  unique() %>% 
  group_by(ids) %>%
  tally()
  
  
corum_annots <- read_tsv("/stor/work/Marcotte/project/rmcox/LECA/ms/gold_stds/corum/allComplexes.txt") %>%
  janitor::clean_names()
