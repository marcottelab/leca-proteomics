library(tidyverse)
library(patchwork)
theme_set(theme_bw(base_size = 12))
pal <- c("#E64B35", "#3578E6", "#4BD2AC", "#E6A435",
         "#665899","#8A8A95", "#353540")

# -----------------------------------
# functions
# -----------------------------------
widget_file_size <- function(p) {
  d <- tempdir()
  withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
  f <- file.path(d, "index.html")
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
}

# -----------------------------------
# main
# -----------------------------------

gene2og <- read_tsv('leca/ppi_ml/data/og_proteomes/nog_mapping/human.euNOG.diamond.mapping.2759') %>% 
  mutate(entry = str_extract(ProteinID,'(?<=\\|)(.*)(?=\\|)')) %>%
  select(entry, ID)

leca_ogs <- read_delim('leca/ppi_ml/annotations/lists/leca_nogs.txt', delim="\n")

df <- read_csv('leca/human_disease/human_disease_groups_manuallyedited.csv') %>%
  select(entry, gene_names_primary, label) %>%
  mutate(label = ifelse(label=="Mental retardation", "Intellectual development\ndisorder", label),
         label = ifelse(label=="Candidiasis", "Familial candidiasis", label)) %>% 
  unique() %>% 
  left_join(gene2og) %>% 
  select(entry, ID, gene_names_primary, label) %>%
  group_by(ID) %>%
  summarize_all(funs(paste(unique(.), collapse = ', '))) %>%
  ungroup() %>% 
  separate_rows(label, sep=", ") %>%
  mutate(leca = ifelse(ID %in% leca_ogs$level_2759, TRUE, FALSE)) %>%
  group_by(label, leca) %>%
  tally() %>%
  arrange(desc(label)) %>%
  pivot_wider(names_from = leca, values_from = n) %>%
  rename(leca_true = `TRUE`, leca_false = `FALSE`) %>%
  replace(is.na(.), as.numeric("0")) %>% 
  mutate(total = leca_true + leca_false) %>%
  filter(total >= 4)

cats <- read_csv('leca/human_disease/disease_categories.021722.csv')

plot <- ggplot(df, aes(x = total, y = (leca_true/total)*100, 
               color = (leca_true/total)*100, label = label)) +
  geom_point(size = 5, alpha = 0.65, stroke = 0.9) +
  scale_color_continuous(type = "viridis", direction = -1,
                         breaks = c(0,50,100)) +
  labs(x = "# of genes associated with each disease\n(total # of OMIM diseases = 295)",
       y = "% of genes tracing back to LECA") +
  theme(legend.title = element_blank(),
        legend.key.height = unit(0.6, "in")) +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5)
plot
prefix = "leca/human_disease/figures/leca_diseases_height5"
h = 5; w = 10.4
plot %>% ggsave(paste0(prefix, ".png"), ., device = "png",
              width = w, height = h, units = "in")
plot %>% ggsave(paste0(prefix, ".pdf"), ., device = "pdf",
              width = w, height = h, units = "in")


# ----- picking labels -----
# library(plotly)
# library(htmlwidgets)
# 
# p <- ggplot(df, aes(x = total, y = leca_true/total, 
#                     color = leca_true/total, label = label)) +
#   geom_point(size = 4, alpha = 0.5, stroke = 0.9) +
#   scale_color_continuous(type = "viridis", direction = -1) +
#   labs(x = "# of genes associated with each disease",
#        y = "% of genes tracing back to LECA") +
#   theme(legend.title = element_blank(),
#         legend.position = "top")
# 
# p_int <- ggplotly(p)
# 
# widget_file_size(p_int)
# saveWidget(p_int, "leca/figures/leca_disease_interactive.html", selfcontained = T, libdir = "lib")
