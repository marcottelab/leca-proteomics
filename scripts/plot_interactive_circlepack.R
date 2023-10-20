library(tidyverse)
library(circlepackeR)
library(data.tree)
library(htmlwidgets)
library(hrbrthemes)
library(treemap)

# ---- STOP AND READ THIS ----
# you have tried many times to get this to work
# you can't, stop spending time here
# learning javascript will be faster
# turn around and go back
# ----------------------------

# functions
widget_file_size <- function(p) {
  d <- tempdir()
  withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
  f <- file.path(d, "index.html")
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
}

# format data
status_col = 'final_status'
workdir <- "/stor/work/Marcotte/project/rmcox/leca/"
clst_annot_file <- paste0(workdir, "ppi_ml/results/walktrap/LinearSVC_100feats_fdr10_4steps_nochloro_dynamic_algo_labels_092523.csv") 

cut_cols = c('cut_676', 'cut_796', 'cut_1202', 'cut_1608', 'cut_2014')

clst_annot <- read_csv(clst_annot_file) %>%
  mutate(human_family_size = str_count(human_gene_names_primary, ',')+1) %>%
  mutate(human_family_size = ifelse(is.na(human_family_size), 0, human_family_size)) %>% 
  mutate(human_genes_fmt = str_replace_all(human_gene_names_primary, pattern = ";", replacement = ",")) %>% 
  mutate(human_genes_fmt = str_replace(human_gene_names_primary, "^([^,]*,[^,]*),.*", "\\1 ...")) %>%
  mutate(human_genes_fmt = str_replace(human_genes_fmt, ", ", "\n")) %>%
  mutate(genes_fmt = ifelse(human_family_size != 0 & human_family_size <= 5, human_genes_fmt, ID)) %>%
  rename(status = status_col) %>% 
  select(ID, human_gene_names_primary, human_genes_fmt, genes_fmt,
         matches("cut_*"), granulated_cmplx_name, status)

pdf <- clst_annot %>%
  filter(!is.na(status), status != "Unclustered") %>% 
  mutate(status_fmt = case_when(status == "Known" ~ "Known protein complex",
                                status == "Novel association" ~ "Novel association",
                                status == "Uncharacterized" ~ "Uncharacterized PPI",
                                TRUE ~ status)) %>%
  select(genes_fmt, any_of(cut_cols), granulated_cmplx_name, status_fmt) %>%
  mutate(root = "leca")

legend_order = c("Known protein complex","Novel association","Uncharacterized PPI",
                 "Large heterogeneous complex")
pdf$status_fmt <- factor(pdf$status_fmt, levels = legend_order)
  

# convert to data.tree structure for circlepackeR
pdf <- pdf %>% 
  mutate(pathString = paste("leca", cut_676, cut_796,
                            granulated_cmplx_name,
                            sep="/"))
leca <- as.Node(pdf)
#p <- circlepackeR(leca, color_min = "hsl(56,80%,80%)", color_max = "hsl(341,30%,40%)")
p <- circlepackeR(leca)
p

outfile <- paste0(workdir, 'figures/interactive_circleplot.html')
widget_file_size(p)
saveWidget(p, outfile, selfcontained = T, libdir = "lib")


# well it works with this url...
url <- "https://gist.githubusercontent.com/mbostock/1093025/raw/05621a578a66fba4d2cbf5a77e2d1bb3a27ac3d4/flare.json"
circlepackeR(url)


# convert to json?
pdf_json <- pdf %>% 
  group_by(root, cut_676, cut_796, cut_1202,
           granulated_cmplx_name) %>%
  summarise(across(everything(), list), .groups = "drop") %>% 
  nest(children = !c(root, cut_676, cut_796,
                           granulated_cmplx_name)) %>% 
  nest(name = !root) %>% 
  toJSON(pretty = TRUE)
pdf_json

write(pdf_json, file=paste0(workdir, 'figures/interactive_circleplot_files/test.json'))

p <- circlepackeR(paste0(workdir, 'figures/interactive_circleplot_files/test.json'))
p

widget_file_size(p)
saveWidget(p, outfile, selfcontained = T, libdir = "lib")