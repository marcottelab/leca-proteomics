# need to concatenate elution vectors from these directories
# need to retain species, experiment, and clade information
dir_animals <- "elut_animals"
dir_plants <- "elut_plants"
dir_tsar <- "elut_tsar"
dir_excavate <- "elut_excavate" 

combine_eluts <- function(data_path, file_pattern, label){
  
  print("Locating data ...")
  files <- dir(data_path, file_pattern) # get file names from given pattern
  
  print("Joining data while retaining filename information ...")
  nested_df <- tibble(filename = files) %>% # create a data frame holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ read_csv(file.path(data_path, .)))) # a new data column
  df <- unnest(nested_df, cols = c(file_contents))
  
  print("Formatting species & experiment info ...")
  df_fmt <- df %>% # extract species and experiment information
    separate(filename, into = c("species", "experiment"), sep="\\.") %>%
    mutate(label = as.character(label))  # add clade label
  
  print("Write out combined results ...")
  write_csv(df_fmt, paste0(data_path, "/", label, str_remove(file_pattern, "\\*")))
  
  print("Done!")
  return(df_fmt)
}

# path to data
workdir <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/"
# pattern to match for files
tidy_norm_pattern <- "*.norm.tidy"

# note: code below takes a long time to run
# files were written out to the data_path with the result
# so can read them back for quicker access
animals <- combine_eluts(data_path = file.path(workdir, dir_animals),
                         file_pattern = tidy_norm_pattern,
                         label = "Amorphea")
plants <- combine_eluts(data_path = file.path(workdir, dir_plants),
                         file_pattern = tidy_norm_pattern,
                         label = "Viridiplantae")
tsar <- combine_eluts(data_path = file.path(workdir, dir_tsar),
                        file_pattern = tidy_norm_pattern,
                        label = "TSAR")
excavate <- combine_eluts(data_path = file.path(workdir, dir_excavate),
                      file_pattern = tidy_norm_pattern,
                      label = "Excavate")
