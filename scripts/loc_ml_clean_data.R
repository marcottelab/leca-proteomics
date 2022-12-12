library(tidyverse)
setwd("/stor/work/Marcotte/project/rmcox/")
palette_pretty <- c("#0072B2","#E69F00","#009E24","#FF0000", "#5530AA")
theme_set(theme_bw(base_size = 12))

human <- read_tsv("leca/localization_ml/data/human_qgo_all.tsv")
yeast <- read_tsv("leca/localization_ml/data/yeast_qgo_all.tsv")
arath <- read_tsv("leca/localization_ml/data/arath_qgo_all.tsv")

clean_go_cc <- function(df, exp){
  
  clean_df <- df %>%
    janitor::clean_names() %>%
    filter(qualifier == "located_in") %>% 
    mutate(Taxon = exp) %>%
    select(gene_product_id:go_name, Taxon) %>%
    unique()
  
    return(clean_df)
    
}

hclean <- clean_go_cc(human, "H. sapiens")
yclean <- clean_go_cc(yeast, "S. cerevisiae")
aclean <- clean_go_cc(arath, "A. thaliana")

all_data <- rbind(hclean, yclean, aclean)

summarized <- all_data %>%
  group_by(go_name, Taxon) %>%
  tally() %>% 
  arrange(desc(Taxon), desc(n))

write_csv(summarized, "leca/localization_ml/results/summarized_cc_counts.csv")

pdata <- summarized %>%
  ungroup %>%
  group_by(Taxon) %>% 
  rename(count = n) %>% 
  slice_max(count, n=13) %>%
  arrange(desc(go_name))

ggplot(pdata, aes(x = go_name, y = count, fill = Taxon)) +
  geom_col() +
  scale_fill_manual(values = palette_pretty) +
  facet_wrap(~Taxon, scales = "free") +
  coord_flip()
  
