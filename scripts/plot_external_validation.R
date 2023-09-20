library(tidyverse)
library(ggrepel)

workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/"
theme_set(cowplot::theme_cowplot())
palette_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85", "#42445A", "#2A303D")
options(scipen=1000)

data <- read_csv(paste0(workdir, 'ext_val_logloss.csv')) %>%
  mutate(clean_labels = case_when(
    set == 'cx' ~ 'HumanNet v3 Co-expression',
    set == 'xlms' ~ 'Bartolec et al. XLMS',
    set == 'y2h' ~ 'Luck et al. Y2H'
    ))
cutoff = 0.437849144630267

ggplot(data, aes(x = mean_ppi_score, y = int_frac*100,
                 color = clean_labels)) +
  geom_line(size = 1) +
  #geom_point(alpha = 0.1) +
  scale_color_manual(values = pal) +
  #scale_y_log10() +
  theme(legend.position = "bottom")


# cumulative proportion plots
data %>%
  group_by(set) %>%
  summarize(int_total = sum(int_frac),
            all_total = sum(total_frac))

ggplot(data, aes(x = mean_ppi_score, y = int_frac, 
                 color = clean_labels)) +
  stat_ecdf(geom="point") +
  labs(y='Cumulative density of overlapping positive pairs',
       x = 'Mean PPI score (bin size = 1000)') +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = cutoff, linetype = "dashed",
             alpha = 0.9) +
  geom_text(label = '10% FDR', x = cutoff+(cutoff*0.25), y = 0.5,
            color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggplot(data, aes(x = mean_ppi_score, y = int_cumu_overlap, 
                 color = clean_labels, group = clean_labels)) +
  geom_point() +
  labs(y='Cumulative fraction of overlapping positive pairs',
       x = 'Mean PPI score (bin size = 1000)') +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = cutoff, linetype = "dashed",
             alpha = 0.75) +
  geom_text(label = '10% FDR', x = cutoff-(cutoff*0.18), y = 1,
            color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = "bottom") -> cumu_plot
cumu_plot

ggsave("leca/ppi_ml/figures/external_ppis_cumulative_overlap.png", cumu_plot, device = "png", width = 8, height = 6, units = "in")

# -----------------------------------------------------------

# density plot
ggplot(data, aes(x = mean_ppi_score, 
                 fill = clean_labels, group = clean_labels,)) +
  geom_density(alpha=0.65) +
  labs(y='Density of overlapping positive pairs',
       x = 'Mean PPI score (bin size = 1000)') +
  scale_fill_manual(values = pal) +
  geom_vline(xintercept = cutoff, linetype = "dashed",
             alpha = 0.75) +
  geom_text(label = '10% FDR', x = cutoff-(cutoff*0.1), y = 7.5,
            color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

# -----------------------------------------------------------

# -----------------------------------------------------------

data <- read_csv(paste0(workdir, 'ext_val_enrichment_10k_bw.csv')) %>%
  mutate(clean_labels = case_when(
    exp == 'cx' ~ 'HumanNet v3 Co-expression',
    exp == 'xlms' ~ 'Bartolec et al. XLMS',
    exp == 'y2h' ~ 'Luck et al. Y2H'
  )) %>%
  mutate(theo_probability = odds_ratio/(odds_ratio+1)) %>%
  mutate(log10_odds = log(odds_ratio, base=10),
         log10_prob = log(prob, base=10))
  

# odds ratio
ggplot(data, aes(x = avg_ppi_score,
                 y = odds_ratio,
                 color = clean_labels)) +
  #geom_point(alpha=0.1) +
  #geom_smooth() +
  geom_line(aes(group = clean_labels), alpha = 0.75,
            size = 2) +
  labs(y='Odds ratio (log scale)',
       x = 'Mean PPI score (bin size = 10,000)') +
  scale_color_manual(values = pal) +
  scale_y_log10() +
  geom_vline(xintercept = cutoff, linetype = "dashed",
             alpha = 0.75) +
  geom_text(label = '10% FDR', x = cutoff-(cutoff*0.1), y = 7.5,
            color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

# probability

calc_r2 <- function(x, y, target){
  
  if(!missing(target)){
  df <- filter(df, exp == target)
  }
  
  fit <- lm(y ~ x)
  r2 <- summary(fit)$adj.r.squared
  pval <- summary(fit)$coef[2,4]
  slope <- fit$coef[[2]]
  intercept <- fit$coef[[1]]
  
  return(r2)
  
}

fit <- mgcv::gam(data$prob ~ data$avg_ppi_score)
summary(fit)

r2_data <- data %>%
  group_by(exp) %>%
  summarize(r2 = calc_r2(avg_ppi_score, prob))

# ------------------------------------------------
# final plot
# ------------------------------------------------
workdir <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/"
theme_set(cowplot::theme_cowplot())

data <- read_csv(paste0(workdir, 'ext_val_enrichment_10k_bw.csv')) %>%
  mutate(clean_labels = case_when(
    exp == 'cx' ~ 'HumanNet v3 Co-expression',
    exp == 'xlms' ~ 'Bartolec et al. XLMS',
    exp == 'y2h' ~ 'Luck et al. Y2H'
  )) %>%
  mutate(theo_probability = odds_ratio/(odds_ratio+1)) %>%
  mutate(log10_odds = log(odds_ratio, base=10),
         log10_prob = log(prob, base=10))

cutoff = 0.437849144630267

pal_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
             "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
             "#7E6148", "#B09C85", "#42445A", "#2A303D")

# reorder legend
data$clean_labels <- factor(data$clean_labels, levels=c('Bartolec et al. XLMS', 
                                                        'Luck et al. Y2H', 
                                                        'HumanNet v3 Co-expression'))

# probability
p <- ggplot(data, aes(x = avg_ppi_score,
                 y = prob,
                 color = clean_labels)) +
  #geom_point(alpha=0.75, size=2) +
  #geom_smooth() +
  geom_line(aes(group = clean_labels), alpha = 0.85, size = 2) +
  # labs(y='Relative likelihood that LECA PPIs are\nrepresented in external data sets (log scale)',
  #      x = 'Mean PPI score (bin size = 10,000)') +
  labs(y='L(external PPIs | LECA PPI score), log scale',
       x = 'Mean PPI score (bin size = 10,000)') +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  scale_color_manual(values = c(pal_npg[2], pal_npg[1], pal_npg[3])) +
  scale_y_log10(limits = c(0.1,100)) +
  #scale_y_continuous(trans="log2") +
  geom_vline(xintercept = cutoff, linetype = "dashed",
             alpha = 0.75) +
  geom_text(label = '10% FDR', x = cutoff-(cutoff*0.13), y = 2.25,
            color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title=element_text(size=12),
        legend.text=element_text(size=10))
p

p %>% ggsave("leca/ppi_ml/figures/likelihood_external_ppis.png", ., device = "png",
                 width = 5, height = 5, units = "in")
p %>% ggsave("leca/ppi_ml/figures/likelihood_external_ppis.pdf", ., device = "pdf",
                 width = 5, height = 5, units = "in")

# ggplot(data, aes(x = prob, y = odds_ratio, color = clean_labels)) +
#   geom_point() +
#   facet_wrap(~clean_labels) +
#   geom_abline()

# -----------------------------------------------------------


# p1 %>% ggsave("leca/ppi_ml/figures/logloss_external_ppis.png", ., device = "png",
#                  width = 10, height = 6, units = "in")

# p %>% ggsave("leca/ppi_ml/figures/logloss_external_ppis.pdf", ., device = "pdf",
#                  width = 6, height = 8, units = "in")

p2 <- data %>%
  group_by(set) %>% 
  ggplot(., aes(x = ppi_score, y = loglike, color = set)) +
  geom_line(size = 2) +
  #geom_point(size = 2) +
  scale_color_manual(values = pal) +
  theme(legend.position = "bottom") +
  #scale_y_log10() +
  facet_wrap(~label)
p2

p2 %>% ggsave("leca/ppi_ml/figures/loglike_external_ppis.png", ., device = "png", width = 10, height = 6, units = "in")
