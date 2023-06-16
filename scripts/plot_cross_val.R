library(tidyverse)
library(patchwork)

dir <- '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cross_val'
files <- list.files(dir, pattern="precision_recall*", full.names = T)

fold_exp = '\\d(?=.csv)'
model_exp = '(?<=recall_).+(?=\\d)'

data <- data_frame(filename = files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest(cols = c(contents)) %>%
  mutate(fold = as.numeric(str_extract(filename, fold_exp))) %>%
  mutate(model = str_extract(filename, model_exp))
data

ap_files <- list.files(dir, pattern="*avg_precision.csv", full.names = T)
ap_data <- data_frame(filename = ap_files) %>%
  mutate(contents = map(filename, ~read_csv(.))) %>%
  unnest() %>%
  mutate(model = str_extract(filename, '(?<=cross_val/).+(?=_GroupKFold)')) %>%
  select(-filename)

ap_stdev <- ap_data %>%
  group_by(model) %>%
  summarize(sd = sd(average_precision))

pdf <- data %>%
  left_join(ap_data, by=c("fold" = "fold", "model" = "model")) %>%
  mutate(ap_label = paste0("(AP=", average_precision, ")")) %>% 
  mutate(legend_label = paste0("GroupKFold ",fold," ",ap_label))

theme_set(theme_bw())
palette_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85")

plot_cv <- function(df, model_name){
  
  textSize = 8
  pointSize = 4
  spaceLegend = 0.1
  
  plot <- df %>%
    filter(model == model_name) %>% 
    group_by(model, fold) %>%
    ggplot(aes(x = recall, y = precision, color = legend_label)) +
    geom_line(size = 1.5) +
    scale_color_manual(values = palette_npg, name = "") +
    theme(legend.position = c(0.35, 0.15),
          legend.title = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          legend.margin = margin(0, 1, 1, 1)) +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    facet_wrap(~model)
    
    return(plot)
  
}

p1 <- plot_cv(pdf, "ExtraTreesClassifier") +
  ylab("Precision")

p2 <- plot_cv(pdf, "LinearSVC") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank())

p3 <- plot_cv(pdf, "SGDClassifier") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank())

patched <- p1 + p2 + p3 & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 0))

# use the tag label as an x-axis label
panel <- wrap_elements(patched) +
  labs(tag = "Recall") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )
panel

panel %>% ggsave("../figures/crossval_x_model_comparison_panel_041723.png", ., device = "png",
                 width = 10, height = 4, units = "in")
panel %>% ggsave("../figures/crossval_x_model_comparison_panel_041723.pdf", ., device = "pdf",
                 width = 10, height = 4, units = "in")
