library(tidyverse)
library(ggrepel)

workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/walktrap/step_sweep"
theme_set(cowplot::theme_cowplot())
pal <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
         "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
         "#7E6148", "#B09C85")

files <- dir(workdir, recursive = T, pattern = "_precision_recall.csv", full.names = T)

fmt_data <- function(x){
  df <- read_csv(x)
  info <- names(df)[1]
  print(info)
  prot_count = str_extract(string = info, pattern = "\\d.+(?=\\))")
  print(prot_count)
  df <- df %>%
    rename(n_clusters = 1) %>%
    mutate(n_prots = prot_count)
  return(df)
}

data <- data_frame(filename = files) %>%
  mutate(contents = map(filename, ~fmt_data(.))) %>%
  unnest(cols = c(contents)) %>%
  mutate(n_steps = str_extract(filename, '(?<=fdr10_).*(?=steps_precision)')) %>%
  mutate(model = str_extract(filename, '(?<=step_sweep/).+?(?=_)')) %>%
  mutate(model_info = paste0(model, "\ntotal # OGs = ", n_prots, "\n10% FDR")) %>% 
  select(-filename)
data

labels <- data %>% 
  group_by(model) %>%
  filter(n_steps == 4 & n_clusters <= 2000) %>%
  unique() %>% 
  slice(seq(1, n(), by = 2))

p <- ggplot(data, aes(x = recall, y = precision, color = model_info)) +
  geom_point(size = 5, alpha = 0.75) +
  scale_color_manual(values = pal) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(),
        legend.title = element_blank()) +
  labs(x = "Recall", y = "Precision") +
  # geom_text_repel(data =  labels, aes(label = n_clusters),
  #                 #point.padding = unit(1, "lines"),
  #                 max.overlaps = Inf,
  #                 box.padding = unit(1.25, "lines"),
  #                 fontface = "bold",
  #                 min.segment.length = unit(0, 'lines'),
  #                 #show.legend = FALSE,
  #                 nudge_y = -0.1
  #                 ) +
  coord_fixed()
p

p %>% ggsave("figures/walktrap_precision-recall.png", ., device = "png",
              width = 8, height = 6, units = "in")
p %>% ggsave("figures/walktrap_precision-recall.pdf", ., device = "pdf",
              width = 8, height = 6, units = "in")
  
library(plotly)
library(htmlwidgets)

plotly_data <- data %>%
  mutate(model_label = paste('model:', model),
         precision_label = paste('precision:', precision),
         recall_label = paste('recall:', recall),
         step_label = paste('# walktrap steps:', n_steps),
         cluster_label = paste('# clusters:', n_clusters),
         hover_label = paste(model_label, precision_label, recall_label,
                             step_label, cluster_label,
                             sep = '\n')) %>%
  filter(n_steps == 4)

p_tmp <- ggplot(plotly_data, aes(x = recall, y = precision, color = model_info, text = hover_label)) +
  geom_point(size = 5, alpha = 0.75) +
  scale_color_manual(values = pal) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin=margin(),
        legend.title = element_blank()) +
  labs(x = "Recall", y = "Precision")

p_int <- ggplotly(p_tmp, tooltip = c("text")) %>%
  layout(
    legend = list(
      orientation = "h", x = 0.1, y = -0.2, 
      title = list(text = '')
    )
  )
p_int

widget_file_size <- function(p) {
  d <- tempdir()
  withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
  f <- file.path(d, "index.html")
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
}

widget_file_size(p_int)

saveWidget(p_int, "figures/walktrap_precision-recall_4steps_interactive.html", selfcontained = T, libdir = "lib")
