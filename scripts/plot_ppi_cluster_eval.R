library(tidyverse)
workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/walktrap/step_sweep"
setwd(workdir)
theme_set(cowplot::theme_cowplot())
pal <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
         "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
         "#7E6148", "#B09C85")

files <- dir(workdir, recursive = T, pattern = "_eval.csv", full.names = T)

fmt_data <- function(file){
  
  df <- read_csv(file) %>%
    mutate(idx = row_number()) %>%
    pivot_longer(!c(cmplx, idx), names_to = "cuts", values_to = "prop") %>% 
    mutate(file = file)
  
  return(df)
  
}

data <- lapply(files, fmt_data) %>%
  bind_rows() %>%
  mutate(n_clst = as.numeric(str_extract(cuts, '(?<=_).*'))) %>% 
  #mutate(feature_set = str_extract(file, '(?<=/).*(?=/)')) %>%
  mutate(model = str_extract(file, '(?<=step_sweep/).+?(?=_)')) %>%
  mutate(walktrap_steps = str_extract(file, '(?<=fdr10_).+?(?=steps_)')) %>%
  mutate(perc = prop*100) %>%
  select(-cuts, -file)

n_obs <- data %>%
  group_by(model, walktrap_steps) %>%
  summarise(count = n_distinct(cmplx))

p1 <- data %>%
  #filter(model == "ExtraTreesClassifier") %>%
  #filter(cmplx < 500) %>%
  #mutate(cmplx = rev(reorder(cmplx, idx))) %>% 
  group_by(cmplx) %>%
    ggplot(aes(x = as.factor(n_clst), y = prop,
               fill = walktrap_steps)) +
      geom_boxplot() +
      #geom_point() +
      scale_fill_brewer(type = "seq", palette = "YlGnBu",
                           name = "# walktrap steps",
                           direction = 1) +
      theme(legend.position="top",
            legend.key.width = unit(0.75, "in"),
            #axis.text.y = element_blank(),
            #axis.ticks.y = element_blank(),
            axis.text.x = element_text(
              angle = 45, vjust = 1, hjust = 1
              )
            ) +
      labs(y = "Fraction of observed gold standard protein complex subunits in the same cluster", 
           x = "# Clusters") +
  facet_wrap(~model, ncol = 1, scales = "free")
p1

p1 %>% ggsave("../figures/walktrap_x_model_comparison_041323.png", ., device = "png",
                 width = 10, height = 10, units = "in")
p1 %>% ggsave("../figures/walktrap_x_model_comparison_041323.pdf", ., device = "pdf",
                 width = 10, height = 10, units = "in")

# summary_stats <- data %>%
#   group_by(model, n_clst, walktrap_steps) %>% 
#   summarize(avg = mean(perc),
#          median = median(perc),
#          stdev = sd(perc),
#          rsd = avg/stdev)
# 
# ggplot(summary_stats, aes(x = as.factor(n_clst),
#                           y = avg)) +
#   geom_col() +
#   theme(axis.text.x = element_text(
#     angle = 45, vjust = 1, hjust = 1
#     )) +
#   facet_wrap(~model, nrow = 1, scales = "free")

files <- dir(workdir, recursive = T, pattern = "steps.csv", full.names = T)

fmt_data <- function(file){
  
  df <- read_csv(file) %>%
    select(ID, matches("*cut*")) %>% 
    pivot_longer(!c(ID), names_to = "cuts", values_to = "cluster") %>%
    mutate(file = file)
  
  return(df)
  
}

data_all <- lapply(files, fmt_data) %>%
  bind_rows() %>%
  mutate(n_clst = as.numeric(str_extract(cuts, '(?<=_).*'))) %>% 
  #mutate(feature_set = str_extract(file, '(?<=/).*(?=/)')) %>%
  mutate(model = str_extract(file, '(?<=step_sweep/).+?(?=_)')) %>%
  mutate(walktrap_steps = str_extract(file, '(?<=fdr10_).+?(?=steps)')) %>% 
  select(-file)

clst_counts <- data_all %>%
  select(-ID, -cuts) %>% 
  group_by(model, walktrap_steps, n_clst, cluster) %>%
  tally()


p2 <- ggplot(clst_counts, aes(x = as.factor(n_clst), y = n,
                        fill = walktrap_steps)) +
  geom_boxplot() +
  scale_fill_brewer(type = "seq", palette = "YlGnBu",
                    name = "# walktrap steps",
                    direction = 1) +
  theme(legend.position="top",
        legend.key.width = unit(0.75, "in"),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.x = element_text(
          angle = 45, vjust = 1, hjust = 1
        )
  ) +
  #scale_y_log10() +
  scale_y_continuous(trans = "log10") +
  labs(y = "Relative protein complex cluster size (# prots in cluster/# clusters)", 
       x = "# Clusters") +
  facet_wrap(~model, ncol = 1, scales = "free")
p2

p2 %>% ggsave("../figures/walktrap_x_model_comparison_cmpx_size_041323.png", ., device = "png",
              width = 10, height = 10, units = "in")
p2 %>% ggsave("../figures/walktrap_x_model_comparison_cmpx_size_041323.pdf", ., device = "pdf",
              width = 10, height = 10, units = "in")

# clst_counts %>%
#   filter(model == "ExtraTreesClassifier") %>%
#   ggplot(aes(x = n,
#              fill = walktrap_steps)) +
#   geom_histogram(bins=100, position="dodge") +
#   scale_fill_brewer(type = "seq", palette = "YlGnBu",
#                     name = "# walktrap steps",
#                     direction = 1) +
#   scale_x_log10()
