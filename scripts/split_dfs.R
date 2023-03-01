input_file <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/input_data/leca_cmplxlist_042922.csv"
cmplx_all <- read_csv(input_file)
cmplx_list <- split(cmplx_all, f = cmplx_all$cmplx_name)
output_dir <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/"

for (i in seq_along(cmplx_list)) {
  filename = gsub(" ", "_", names(cmplx_list)[i])
  filename_fmt = paste0(output_dir, filename, ".csv")
  write_csv(cmplx_list[[i]], filename_fmt)
}
