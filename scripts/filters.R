library(tidyverse)
library(data.table)
theme_set(cowplot::theme_cowplot())
options(scipen=10000)

workdir <- "/stor/work/Marcotte/project/rmcox/LECA/ms"
palette_pretty <- c("#009E24", "#0072B2","#E69F00", "#FF0000",
                    "#979797", "#5530AA", "#1E1E1E")

# read in data
exp_info <- read_csv(file.path(workdir, "cfms/meta/leca_experiments.csv")) %>%
  janitor::clean_names()
clades <- read_csv(file.path(workdir, "meta/speciesinfo_clades.csv")) %>%
  select(code, clade)
totals <- read_tsv(file.path(workdir, "cfms2/results/leca_euk_elut.dollonogs.totalcounts.txt"))
elut_matrix <- read_tsv(file.path(workdir, "cfms2/results/leca_euks_elut.filtdollo.raw.txt"))
fmat_dt <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_firstpass/model_training/featmat_labeled"),
                             sep = ",")
obs_prots <- elut_matrix %>%
  select(1)

# -------------------------------------------
# choose filter for total PSMs across all experiments
# -------------------------------------------

fmat_dt_fmt <- fmat_dt[, c("ID1", "ID2") := data.table::tstrsplit(ID, " ", fixed=TRUE)][]
fmat_labels <- fmat_dt_fmt %>%
  select(ID1, ID2, label) %>%
  mutate(interaction = case_when(label == 1 ~ "Positive",
                                 label == 0 ~ "Other",
                                 label == -1 ~ "Negative"))

fmat_labels_totals <- fmat_labels %>%
  left_join(totals, by = c("ID1" = "X1")) %>%
  rename(ID1_PSMs = total_PSMs) %>%
  left_join(totals, by = c("ID2" = "X1")) %>%
  rename(ID2_PSMs = total_PSMs) %>%
  filter(!is.na(ID1_PSMs)) %>%
  filter(!is.na(ID2_PSMs)) %>%
  filter(interaction %in% c("Positive", "Negative"))

cutoff = 150
psm_cf <- ggplot(fmat_labels_totals, aes(x = ID2_PSMs, y = ID1_PSMs,
                               color = interaction)) +
  geom_point(alpha = 0.3) +
  geom_vline(xintercept = cutoff) +
  geom_hline(yintercept = cutoff) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c(palette_pretty[3],
                                palette_pretty[2])) +
  facet_wrap(~interaction)
psm_cf

ggsave(file.path(workdir, "../figures/PSM_cutoff.png"), psm_cf,
       device = "png", width = 8, height = 4, units = "in")
ggsave(file.path(workdir, "../figures/PSM_cutoff.pdf"), psm_cf, 
       device = "pdf", width = 8, height = 4, units = "in")


# -------------------------------------------
# make clade-specific filtered featmats
# -------------------------------------------

# function for writing files that divide data into clades
write_clade <- function(df, out){
  
  new_df <- df %>%
    select(code) %>%
    mutate(code = (tolower(code))) %>%
    write_delim(., out, col_names = F, delim = "\n")
  return(new_df)
  
}

# function for concatenating elut data
concat_elut <- function(elut_list, psm_thres, outfile_prefix){
  
  # function for reading in data
  read_elut <- function(elut_file){
    
    elut_df <- read_delim(elut_file, delim = ",")
    
    return(elut_df)
  }
  
  results <- mapply(read_elut, 
                    elut_file = elut_list)
  results[1] <- NULL
  
  # takes a long time
  elut_concat <- results %>%
    reduce(full_join, by = 'X1')
  elut_concat[is.na(elut_concat)] <- 0
  
  path <- dirname(normalizePath(elut_list[1]))
  
  # these take a long time to write
  write_tsv(elut_concat, paste0(path, "/", outfile_prefix, ".raw.txt"))
}

# function for reading in cfms data, with options for filtering and normalization
parse_elut <- function(elutfile, filtfile = NULL, norm = FALSE, 
                       write = FALSE, outfile_suffix = ".parse.elut",
                       sep = "\t"){
  
  elut_df <- read_delim(elutfile, delim = sep)
  
  if(!missing(filtfile)){
    
    filter <- read_delim(filtfile, delim = "\n", col_names = FALSE)
    elut_df <- elut_df %>%
      filter(X1 %in% filter$X1)
  }
  
  if(norm == TRUE){
    
    elut_matrix <- elut_df %>% 
      select(-X1)
    
    elut_norm <- 
      t(apply(elut_matrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
    
    elut_out <- data.frame(elut_norm, orthogroup = elut_df$X1)
    
    elut_df <- elut_out %>%
      select(orthogroup, everything())
  }
  
  if(write == TRUE){
    
    filename <- tools::file_path_sans_ext(elutfile)
    outfilename <- paste0(filename, outfile_suffix)
    write_csv(elut_df, outfilename)
    
  } 
  
  return(elut_df)
  
}

animals <- filter(clades, clade == "Amorphea")
plants <- filter(clades, clade == "Archaeplastida")
tsar <- filter(clades, clade == "TSAR")
excavate <- filter(clades, clade == "Excavate")

# specify which species belong to each clades
write_clade(animals, file.path(workdir, "cfms2/animal_codes.txt"))
write_clade(plants, file.path(workdir, "cfms2/plant_codes.txt"))
write_clade(tsar, file.path(workdir, "cfms2/tsar_codes.txt"))
write_clade(excavate, file.path(workdir, "cfms2/excavate_codes.txt"))

# on the CLI: used codes to sort different clades into different directories:
# "cfms2/elut_animals"
# "cfms2/elut_plants"
# "cfms2/elut_tsar"
# "cfms2/elut_excavate"

# locate data for each clade, filter, concatenate, and write out
# animals
elut_animals = file.path(workdir, "cfms2/elut_animals")
elutlist_animals <- dir(elut_animals, pattern = "*filtdollo.elut",
                  full.names = TRUE)
concat_elut(elut_list = elutlist_animals, 
            outfile_prefix = "animals_concat")
parse_elut(elutfile = file.path(workdir, "cfms2/elut_animals/animals_concat.raw.txt"),
           filtfile = file.path(workdir, "cfms2/filt150p_ids.txt"),
           write = TRUE,
           outfile_suffix = ".150p.elut")
parse_elut(elutfile = file.path(workdir, "cfms2/elut_animals/animals_concat.raw.150p.elut"),
           filtfile = file.path(workdir, "gold_stds/ppis/all.euNOG.gold_prots.unique.txt"),
           write = TRUE,
           sep = ",",
           outfile_suffix = ".gold.elut")

# plants
elut_plants = file.path(workdir, "cfms2/elut_plants")
elutlist_plants <- dir(elut_plants, pattern = "*filtdollo.elut",
                        full.names = TRUE)
concat_elut(elut_list = elutlist_plants,
           outfile_prefix = "plants_concat")
parse_elut(elutfile = file.path(workdir, "cfms2/elut_plants/plants_concat.raw.txt"),
           filtfile = file.path(workdir, "cfms2/filt150p_ids.txt"),
           write = TRUE,
           outfile_suffix = ".150p.elut")
parse_elut(elutfile = file.path(workdir, "cfms2/elut_plants/plants_concat.raw.150p.elut"),
           filtfile = file.path(workdir, "gold_stds/ppis/all.euNOG.gold_prots.unique.txt"),
           write = TRUE,
           sep = ",",
           outfile_suffix = ".gold.elut")


# tsar
elut_tsar = file.path(workdir, "cfms2/elut_tsar")
elutlist_tsar <- dir(elut_tsar, pattern = "*filtdollo.elut",
                       full.names = TRUE)
concat_elut(elut_list = elutlist_tsar,
           outfile_prefix = "tsar_concat")
# only IDs sampled at 150 PSMs
parse_elut(elutfile = file.path(workdir, "cfms2/elut_tsar/tsar_concat.raw.txt"),
           filtfile = file.path(workdir, "cfms2/filt150p_ids.txt"),
           write = TRUE,
           outfile_suffix = ".150p.elut")
# only gold standard IDs
parse_elut(elutfile = file.path(workdir, "cfms2/elut_tsar/tsar_concat.raw.150p.elut"),
           filtfile = file.path(workdir, "gold_stds/ppis/all.euNOG.gold_prots.unique.txt"),
           write = TRUE,
           sep = ",",
           outfile_suffix = ".gold.elut")

# excavate
elut_excavate = file.path(workdir, "cfms2/elut_excavate")
elutlist_excavate <- dir(elut_excavate, pattern = "*filtdollo.elut",
                       full.names = TRUE)
concat_elut(elut_list = elutlist_excavate,
           outfile_prefix = "excavate_concat")
parse_elut(elutfile = file.path(workdir, "cfms2/elut_excavate/excavate_concat.raw.txt"),
           filtfile = file.path(workdir, "cfms2/filt150p_ids.txt"),
           write = TRUE,
           outfile_suffix = ".150p.elut")
parse_elut(elutfile = file.path(workdir, "cfms2/elut_excavate/excavate_concat.raw.150p.elut"),
           filtfile = file.path(workdir, "gold_stds/ppis/all.euNOG.gold_prots.unique.txt"),
           write = TRUE,
           sep = ",",
           outfile_suffix = ".gold.elut")


# -------------------------------------------
# filter by minimum correlation coefficient
# -------------------------------------------

# DONE --- applied the filter using python
# python3 /stor/work/Marcotte/project/rmcox/LECA/scripts/filter_clades.py -f ids_exps_species.raw.150p.p3.csv -t 2
# these are the results
# p3c2_res <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/features/ids_exps_species.raw.150p.p3.clade_filt.csv"),
#                               sep = ",", col.names = T)
# p3c2_res <- read_csv(file.path(workdir, "cfms2/cfmsflow_012022/features/ids_exps_species.raw.150p.p3.clade_filt.csv"))
# p3c2_res <- p3c2_res %>%
#   select(ID) %>%
#   unique()
# write_delim(p3c2_res, file.path(workdir, "cfms2/cfmsflow_012022/feature_filters/filtp3c2_ids.txt"), col_names = F, delim = "\n")

# function for filtering feature data
parse_feats <- function(featfile, filtfile = NULL,
                        write = FALSE, outfile_suffix = ".parse.feat",
                        sep = ","){
  
  new_col1 <- str_extract(basename(featfile), pattern = "^(.*?\\..*?)\\.")
  new_col2 <- str_extract(basename(featfile), pattern = "([^\\.]+)(?:\\.[^\\.]+?){2}$")
  new_col_final <- paste0(new_col1, new_col2)
  
  feat_df <- data.table::fread(featfile, sep = sep, header = T) %>%
    rename_at(vars(-ID1, -ID2), ~paste0(new_col1, new_col2))
  feat_df[,ID:=paste(ID1,ID2,sep=" ")]
  feat_df[,ID1:=NULL]
  feat_df[,ID2:=NULL]
  
  feat_df <- select(feat_df, ID, everything())
  print(feat_df)
  
  if(!missing(filtfile)){
    
    filter <- read_delim(filtfile, delim = "\n", col_names = FALSE)
    feat_df <- feat_df %>%
      filter(ID %in% filter$X1)
  }
  
  if(write == TRUE){
    
    filename <- tools::file_path_sans_ext(featfile)
    outfilename <- paste0(filename, outfile_suffix)
    write_csv(feat_df, outfilename)
    
  } 
  
  return(feat_df)
  
}


# DONE --- parsed files written out to file.path(workdir, "cfms2/cfmsflow_012022/features/")
filtp3c2 <- read_csv(file.path(workdir, "cfms2/cfmsflow_012022/feature_filters/filtp3c2_ids.txt"),
                     col_names = F)
data_path <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/features"   # path to the data
files <- dir(data_path, pattern = "*.feat$", full.names = T) # get file names

system.time(
p3c2_parse <- mapply(parse_feats,
                     featfile = files,
                     filtfile = file.path(workdir, "cfms2/cfmsflow_012022/feature_filters/filtp3c2_ids.txt"),
                     write = TRUE,
                     outfile_suffix = ".p3c2.feat")
)

system.time({
  big_df = foreach(i = fn_test,
                   .packages = c("data.table", "tidyverse")) %dopar% {
                     fread(i, header = TRUE)
                   } %>%
    reduce(full_join, by = "ID")
})


p3c2_parse[1] <- NULL
p3c2_concat <- p3c2_parse %>%
  reduce(full_join, by = 'ID')
p3c2_concat[is.na(p3c2_concat)] <- 0


########## let's try something weird

library(doParallel)
#p_load(doParallel,data.table,stringr)

# get the file name
data_path <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/features"
dir(data_path, pattern = "*p3c3*", full.names = T) -> fn
fn_test <- sample(fn, 10)

# use parallel setting
(cl = detectCores() %>% 
    makeCluster()) %>% 
  registerDoParallel()

# read and bind
system.time({
  big_df = foreach(i = fn,
                   .packages = c("data.table", "tidyverse")) %dopar% {
                     fread(i, header = TRUE)
                   } %>%
    reduce(full_join, by = "ID")
})

# end of parallel work
parallel::stopCluster(cl)
##########

big_df[is.na(big_df)] <- 0
system.time(
  data.table::fwrite(big_df, 
                     file.path(workdir, 
                               "cfms2/cfmsflow_012022/featmat/allexps_p3c2_featmat"), 
                     sep = ",")
)
print("Done!")


# -------------------------------------------
# final feature matrix
# -------------------------------------------
# filter for proteins with > 150 total PSMs
filt150p <- read_tsv(file.path(workdir, "cfms2/filt150p_ids.txt"), col_names = F)
goldstds <- read_delim(file.path(workdir, "gold_stds/all.gold.euNOGs.txt"),
                      delim = "\n",
                      col_names = F)

# filtered feature matrix (dollo, gold std, 0.3 pearson's R + 2 species, 150 PSMs)
fmat_pRfilt_dt <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_120721/featmat/pRfilt_p3_s2.csv"),
                             sep = ",")
fmat_pRfilt_dt[, c("ID1", "ID2") := data.table::tstrsplit(ID, " ", fixed=TRUE)][]
fmat_pR_150p <- fmat_pRfilt_dt %>%
  filter(ID1 %in% filt150p$X1 & ID2 %in% filt150p$X1)  # PSM filter
fmat_pR_150p[,ID:=paste(ID1,ID2,sep=" ")]
fmat_pR_150p[,ID1:=NULL]
fmat_pR_150p[,ID2:=NULL]

# selected external humap features
fmat_humap <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_120721/featmat/humap_ext_feats_final.csv"),
                                sep = ",")
fmat_humap[,ID:=paste(acc1,acc2,sep=" ")]
fmat_humap[,acc1:=NULL]
fmat_humap[,acc2:=NULL]
fmat_humap_grouped <- fmat_humap %>%
  group_by(ID) %>% 
  summarize_all(max)
fmat_humap_grouped <- as.data.table(fmat_humap_grouped)
data.table::fwrite(fmat_humap_grouped, 
                   file.path(workdir, "cfms2/cfmsflow_120721/featmat/featmat_humap"),
                  sep = ",")

# feature matrices for whole data set, and by clade
fmat_allconcat <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_allconcat"),
                                  sep = ",")
fmat_animals <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_animals"),
                                  sep = ",")
fmat_plants <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_plants"),
                                 sep = ",")
fmat_tsar <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_tsar"),
                               sep = ",")
fmat_excavate <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_excavate"),
                                   sep = ",")
fmat_humap <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_humap"),
                                sep = ",")
fmat_p3c2 <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_012022/featmat/featmat_allexps_p3c2"),
                               sep = ",")

joined_feats <- fmat_humap[fmat_p3c2, on = 'ID']
joined_feats <- fmat_allconcat[joined_feats, on = 'ID']
joined_feats <- fmat_animals[joined_feats, on = 'ID']
joined_feats <- fmat_plants[joined_feats, on = 'ID']
joined_feats <- fmat_tsar[joined_feats, on = 'ID']
joined_feats <- fmat_excavate[joined_feats, on = 'ID']

joined_feats[is.na(joined_feats)] <- 0

system.time(
  data.table::fwrite(joined_feats, 
                     file.path(workdir, 
                               "cfms2/cfmsflow_012022/featmat/featmat_final_lj"), 
                     sep = ",")
)


# -------------------------------------------
# results
# -------------------------------------------
res <- read_csv(file.path(workdir, "cfms2/cfmsflow_120721/results/extratrees/clustering_fdr25.csv"))

res_goldstd <- filter(res, ID %in% goldstds$X1)

# all result proteins are gold std proteins...
# were all input proteins gold std?
humap_ids <- data.table::fread(file.path(workdir, "cfms2/cfmsflow_120721/featmat/humap_ext_feats_final.csv"),
                                sep = ",") %>%
  select(acc1, acc2)

humap_ids <- data.frame(ids=unlist(humap_ids, use.names = FALSE))
humap_ids <- unique(humap_ids)

# nope there are way more humap IDs than gold standard IDs

# I can do this a little more directly with the scored interactions file
scored <- read_tsv(file.path(workdir, "cfms2/cfmsflow_120721/model_training/scored_interactions_extratrees")) %>%
  select(ID) %>%
  separate(ID, into = c("ID1", "ID2"), sep = " ")

scored_ids <- data.frame(ids=unlist(scored, use.names = FALSE))
scored_ids <- unique(scored_ids)

scored_ids_gold <- scored_ids %>%
  filter(ids %in% goldstds$X1)

# -------------------------------------------
# final feature matrix based on top features from RFE model ?
# -------------------------------------------
topfeats <- read_delim(file.path(workdir, "cfms2/cfmsflow_120721/model_training/topfeats_rfe.txt"),
                                 delim = "\n")

top_humap <- fmat_humap_grouped %>%
  select(ID, any_of(topfeats$feature))
         