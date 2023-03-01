trim_wspace <- function(x) gsub("^\\s+|\\s+$", "", x)

fmt_subcc <- function(clst_df){
  
  df_subcc_fmt <- clst_df %>%
    select(1:25, human_subcellular_location_cc) %>%
    mutate(localization = str_extract_all(human_subcellular_location_cc,
                                          '(?<=: )(.+?)(?=.$| \\{ECO)')) %>%
    unnest(cols = c(localization)) %>% 
    mutate(localization = str_remove(localization, '\\[.*\\]:')) %>%
    mutate(localization = str_remove(localization, ';.*')) %>%
    mutate(localization = str_remove(localization, 'Note=.*')) %>%
    separate_rows(localization, sep = '\\. ') %>%
    mutate(localization = str_remove(localization, ',\\s.*')) %>%
    mutate(localization = str_remove(localization, '\\.$')) %>%
    mutate(localization = trim_wspace(localization)) %>%
    mutate(localization = case_when(startsWith(as.character(human_gene_names_primary), "RPL") ~ "Multi-compartment",
                                    startsWith(as.character(human_gene_names_primary), "RPS") ~ "Multi-compartment",
                                    startsWith(as.character(human_gene_names_primary), "EIF2B") ~ "Cytoplasm",
                                    TRUE ~ localization)) %>% 
    mutate(localization = case_when(str_detect(localization, 'Mitochondrion') == TRUE ~ 'Mitochondrion',
                                    str_detect(localization, '[rR]eticulum|Microsome') == TRUE ~ 'Endoplasmic reticulum',
                                    str_detect(localization, 'Lysosome') == TRUE ~ 'Lysosome',
                                    str_detect(localization, 'Cytoplasm$') == TRUE ~ 'Cytoplasm',
                                    str_detect(localization, '[cC]ytosol$') == TRUE ~ 'Cytoplasm',
                                    str_detect(localization, 'Nucleus|Chromosome') == TRUE ~ 'Nucleus',
                                    str_detect(localization, '[nN]ucleolus') == TRUE ~ 'Nucleus',
                                    str_detect(localization, 'Early endosome|Endosome membrane|Late endosome|Endosome') == TRUE ~ 'Endomembrane',
                                    str_detect(localization, '[mM]embrane|receptor|surface') == TRUE ~ 'Membrane',
                                    str_detect(localization, 'Secreted') == TRUE ~ 'Secreted',
                                    str_detect(localization, 'autophago|granule') == TRUE ~ 'Cytoplasm',
                                    TRUE ~ localization)) %>%
    mutate(secreted = case_when(localization == "Secreted" ~ "TRUE",
                                TRUE ~ "FALSE")) %>% 
    group_by(group_by(across(c(-localization)))) %>%
    summarize(loc_grouped = list(unique(localization))) %>%
    mutate(num_compartments = lengths(loc_grouped)) %>%
    mutate(loc_fmt = case_when(is.na(loc_grouped) ~ "Unknown",
                               loc_length > 1 ~ "Multi-compartment",
                               TRUE ~ as.character(loc_grouped))) %>%
    mutate(loc_fmt = case_when(loc_fmt == "Cell projection" ~ "Membrane",
                               loc_fmt == "Peroxisome" ~ "Endomembrane",
                               loc_fmt == "Lysosome" ~ "Endomembrane",
                               loc_fmt == "Recycling endosome" ~ "Endomembrane",
                               loc_fmt == "Cytoplasmic vesicle" ~ "Endomembrane",
                               loc_fmt == "Dynein axonemal particle" ~ "Cytoskeleton",
                               loc_fmt == "Cell projection" ~ "Cytoskeleton",
                               loc_fmt == "Cell junction" ~ "Cytoskeleton",
                               loc_fmt == "Midbody" ~ "Cytoskeleton",
                               loc_fmt == "Golgi apparatus" ~ "Endomembrane",
                               loc_fmt == "Endoplasmic reticulum" ~ "Endomembrane",
                               loc_fmt == "Secreted" ~ "Membrane",
                               TRUE ~ loc_fmt)) %>% 
    ungroup()
  
  df_subcc_fmt <- df_subcc_fmt %>%
    select(ID, loc_fmt, num_compartments, secreted)
  
  left_join(clst_df, df_subcc_fmt)
  
  return(df_fmt)
}

fmt_genes <- function(clst_df){
  
  df_fmt <- clst_df %>%
    mutate(gene_name = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ..."))
  
}

df <- read_csv("/project/rmcox/LECA/results/clustering_fdr15_noribo.fmt.csv") %>%
  select(1:27)

df_fmt <- fmt_subcc(df)

leca_clst_fmt <- df_fmt %>%
  mutate(gene_name = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>% 
  mutate(phylo_sums = rowSums(lclst_filt[11:14]))
# 
# write_tsv(leca_clst_fmt, "/project/rmcox/LECA/results/leca_clst_annot.fmt.03292022.tsv")