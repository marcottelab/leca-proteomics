library(tidyverse)

# scores
ppi_scores <- read_tsv("/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/scored_interactions_top50k") 

# annotations
annots <- read_tsv("/stor/work/Marcotte/project/rmcox/LECA/ms/annotations/leca_euNOGs_annotated.tsv") %>%
  select(ID, human_gene_names_primary, human_protein_names,
         arath_gene_names_primary, arath_protein_names)
annot_cols <- c("human_gene_names_primary", "human_protein_names", "arath_gene_names_primary", "arath_protein_names")

# test train status
testfile <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.test_ppis.ordered"
trainfile <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered"

test <- readr::read_csv(testfile, col_names = F) %>%
  mutate(set = "test") %>%
  rename(ID = X1)
train <- readr::read_csv(trainfile, col_names = F) %>%
  mutate(set = "train") %>%
  rename(ID = X1)

test_train <- rbind(test, train)

ppi_scores_set <- left_join(ppi_scores, test_train)
ppi_scores_set[is.na(ppi_scores_set)] <- "predicted"

ppi_scores_annotated <- ppi_scores_set %>%
  separate(ID, into = c("ID1", "ID2"), sep = " ") %>% 
  left_join(annots, by = c("ID1" = "ID")) %>%
  rename_with(.cols = all_of(annot_cols), function(x){paste0("ID1_", x)}) %>%
  left_join(annots, by = c("ID2" = "ID")) %>%
  rename_with(.cols = all_of(annot_cols), function(x){paste0("ID2_", x)})

ppi_scores_annotated[] <- lapply(ppi_scores_annotated, str_trunc, 25000, ellipsis = "...")

write_csv(ppi_scores_annotated,
          "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/results_lj.noribo/annotatedscores_top50k.fmt.csv")  
