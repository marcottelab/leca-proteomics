library(tidyverse)
library(textshape)

hclust_elut <- function(concat_elut_file, outfile, normalize == TRUE){
  
  elut_concat <- read_tsv(concat_elut_file)
  
  # normalization
  if(normalize){
  
  norm01 <- function(x){(x-min(x))/(max(x)-min(x))}
  print("normalization in progress...")
    
  norm <- elut_concat %>%
    select(-X1, -total_PSMs) %>%
    norm01()
  
  elut_norm <- data.frame(norm, 
                          orthogroup = elut_concat$X1) %>%
    tibble::column_to_rownames(var="orthogroup")
  print("normalization complete.")
  
  # cluster
  print("hierarchical clustering in progress...")
  elut_cluster <- elut_norm %>%
    cluster_matrix(dim = "row")
  print("clustering complete.")
  
  } else {
  
  print("normalization set to FALSE, clustering raw matrix...")
  
    # cluster
  elut_cluster <- elut_concat %>%
    select(-total_PSMs) %>% 
    tibble::column_to_rownames(var = "X1") %>%
    cluster_matrix(dim = "row")
  
  print("clustering complete.")
  
  }
  
  print("writing results to", outfile)
  write.table(elut_cluster, outfile)
  
}

#hclust_elut(concat_elut_file, outfile)
  
  