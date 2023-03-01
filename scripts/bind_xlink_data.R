xlink <- read_csv("/project/dy4652/cilia_xlink/result.csv") %>%
  select(X1, ID1, ID2, is_crosslinked, P_1, set, matches("ID1_human*"), matches("ID2_human*")) %>%
  filter(is_crosslinked == TRUE)

write_csv(xlink, "/project/rmcox/LECA/results/leca-ppis_tett-cilia_xlink.csv")

data_path <- "/project/rmcox/LECA/xlink/raw_data/"   # path to the data
files <- dir(data_path, pattern = "*.csv") # get file names
all_xlink_data <- data_frame(filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(data_path, .))) # a new data column
  ) %>%
  unnest()

all_xlink_data

interxlink <- filter(all_xlink_data, Protein1 != Protein2) %>%
  select(Protein1, Protein2) %>%
  unique()

write_csv(all_xlink_data, "/project/rmcox/LECA/xlink/processed_data/all_xlinks_032222.csv")

