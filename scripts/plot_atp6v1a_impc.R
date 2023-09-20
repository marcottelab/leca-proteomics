library(tidyverse)
library(ggsignif)
theme_set(theme_bw(base_size = 12))
pal <- c("#56B4E9", "#F0E442", "#009E24", "#E69F00", 
         "#FF0000", "#979797", "#5530AA", "#1E1E1E")
pal_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85")

data_file = 'leca/human_disease/data/Atp6v1a_IMPC_1201780.tsv'
mutant <- expression("Atp6v1a"^paste(em1(IMPC)/"+"))
#expression(paste("Atp6v1a"^"em1(IMPC)/+", " EP")))

df <- read_tsv(data_file) #%>%
  #mutate(label = ifelse(biological_sample_group == "control", "Control", paste0("~Atp6v1a^{em1(IMPC)/+}")))
  #mutate(pval = ifelse(biological_sample_group == "female", "p=3.33×10-05", ))


p <- ggplot(df, aes(y = data_point, x = sex)) +
  #geom_violin(alpha = 0.8) +
  geom_boxplot(alpha = 0.9, aes(fill = biological_sample_group, color = biological_sample_group)) +
  #geom_jitter() +
  scale_fill_manual(values = c(pal_npg[5], pal_npg[6]),
                    labels = c("Control", mutant)) +
  scale_color_manual(values = c(pal_npg[1], pal_npg[4]),
                     labels = c("Control", mutant)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank()) +
  labs(y = "Bone Mineral Content (g)") +
  geom_signif(stat="identity", parse=TRUE,
              data=data.frame(x=c(0.8, 1.85), xend=c(1.2, 2.2),
                              y=c(0.85, 0.8), 
                              annotation = c("p=3.33×10^-5", "p=0.00813")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation),
              textsize=2, tip_length=0) +
  ylim(c(0.3,0.9))

p %>% ggsave("leca/human_disease/figures/atp6v1a_impc.png", ., device = "png", 
             width = 3.5, height = 3.5, units = "in")
p %>% ggsave("leca/human_disease/figures/atp6v1a_impc.pdf", ., device = "pdf", 
             width = 3.5, height = 3.5, units = "in")
