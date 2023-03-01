

file_path <- "/project/rmcox/LECA/ms/gold_stds/ppis/"
file_pattern <- "*annotations*"
annot_files <- dir(file_path, pattern = file_pattern, full.names = TRUE)
annot_files

annot_data <- annot_files %>%
  map(read_csv) %>%    # read in all the files individually
  reduce(rbind)        # reduce with rbind into one dataframe
data

tax_ids <- read_csv("/project/rmcox/LECA/ms/gold_stds/gold_std_info.csv") %>%
  select(tax_code, species, species_long, taxa)

annot_fmt <- annot_data %>%
  left_join(tax_ids, by = c("taxonomy_identifier" = "tax_code")) %>%
  select(ids_fmt, number_complex_ac, recommended_name,
         aliases_for_complex, taxonomy_identifier, species,
         species_long, taxa, everything())

write_csv(annot_fmt, paste0(file_path, "all.annotations.gold.cmplx.csv"))
