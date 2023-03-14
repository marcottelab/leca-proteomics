path = "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign_totals/"
files = list.files(path, pattern="*csv")

fmt_data <- function(file) {
  
  species = str_split_i(file, "_", 1)
  df <- read_csv(file, col_names = T) %>%
    select(protein, peptide, matches('*_total')) %>%
    group_by(protein) %>%
    summarize_if(is.numeric, sum, na.rm=TRUE)
  
  return(df)
  
}

df <- list.files(path="/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign_totals/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows

df = 


