library(tidyverse)

label_featmat <- function(feat_file, ppis_file, out_file){
  
  print("Reading in feature matrix ... ")
  featmat_dt <- data.table::fread(feat_file, sep = ",")  # read in feature matrix
  featmat_dt[is.na(featmat_dt)] <- 0  # convert missing values to 0
  
  print("Reading in gold standard protein assemblies ... ")
  # read in gold standard complexes
  gold_ppis <- readr::read_delim(ppis_file, delim = "\n", col_names = FALSE) %>%
    dplyr::mutate(complex = row_number()) %>%
    dplyr::mutate(member_id = strsplit(as.character(X1), " ")) %>% 
    tidyr::unnest(member_id) %>%
    dplyr::group_by(member_id)  %>%
    dplyr::summarise(allplex = list(as.numeric(unique(complex))))
  
  gold_ppis <- data.table::as.data.table(gold_ppis)
  
  # set up comparison columns for labeling
  featmat_dt[, c("ID1", "ID2") := data.table::tstrsplit(ID, " ", fixed=TRUE)][]
  
  compare_df <- as.data.frame(featmat_dt[, c('ID1', 'ID2')])
  
  print("Evaluating true/false interactions ... ")
  # compute labels (1 = positive, -1 = negative)
  system.time(
    labeled_df <- compare_df %>%
      rowwise() %>% 
      mutate(id1_loc = match(ID1, gold_ppis$member_id)) %>%
      mutate(id1_vec = gold_ppis[id1_loc, 2]) %>%
      mutate(id2_loc = match(ID2, gold_ppis$member_id)) %>%
      mutate(id2_vec = gold_ppis[id2_loc, 2]) %>%
      mutate(eval = list(intersect(unlist(id1_vec), unlist(id2_vec)))) %>%
      mutate(eval_len = length(eval)) %>% 
      mutate(label = case_when(eval_len > 0 ~ 1,
                               TRUE ~ -1))
    
  )
  print("Evaluation matrix:")
  print(labeled_df)
  
  print("Formatting labeled matrix ...")
  # extract label column
  labels <- select(labeled_df, label)
  
  # tack it onto the featmat
  system.time( 
    final_dt <- cbind(featmat_dt, labels)
  )
  
  print("Writing final matrix ...")
  # write it out
  system.time(
  data.table::fwrite(final_dt, out_file, sep = ",")
  )
  print("Done!")
}

feat_file <- "/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_norm/featmat"
ppis_file <- "/project/rmcox/LECA/ms/gold_stds/all.gold.cmplx.merged.humap.txt"
out_file <- "/project/rmcox/LECA/ms/cfms2/cfmsflow/feat_selection/featmat.gold.norm.labeled.csv"

label_featmat(feat_file, ppis_file, out_file)
