library(tidyverse)

workdir <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/"
setwd(workdir)

pr1f <- "cfmsflow_firstpass/model_training/precisionrecall_extratrees"
pr2f <- "cfmsflow_012022/results_lj/precisionrecall_fdr15"
pr3f <- "cfmsflow_012022/results_lj.noribo/precisionrecall_fdr15_noribo"


pr1 <- readr::read_csv(pr1f) %>%
  filter(label == "test") %>%
  mutate(set = "Iteration #1") %>%
  select(Recall, Precision, set)

pr2 <- readr::read_csv(pr2f) %>%
  filter(label == "test") %>%
  mutate(set = "Iteration #2") %>%
  select(Recall, Precision, set)

pr3 <- readr::read_csv(pr3f) %>%
  filter(label == "test") %>%
  mutate(set = "Iteration #3") %>%
  select(Recall, Precision, set)

all_pr <- rbind(pr1, pr2, pr3)

ggplot(all_pr, aes(x = Recall, y = Precision, color = set, group = set)) +
  geom_line()
