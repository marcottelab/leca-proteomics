
cilia <- read_csv("/stor/work/Marcotte/project/rmcox/leca/xlink/david/cilia_xlink_new/results/annotatedscores_top50k_isCrosslinked2_all_xlinks.csv") %>%
  rename(cilia_xlink = is_crosslinked)

bodies <- read_csv("/stor/work/Marcotte/project/rmcox/leca/xlink/david/cilia_xlink_new/results/annotatedscores_top50k_isCrosslinked2_BodiesXL_1-3a-g.csv") %>%
  rename(bodies_xlink = is_crosslinked) %>%
  select(ID1, ID2, bodies_xlink)

all_xlink <- left_join(cilia, bodies, by=c("ID1", "ID2")) %>%
  select(ID1, ID2, cilia_xlink, bodies_xlink, everything()) %>%
  arrange(desc(bodies_xlink), desc(cilia_xlink), desc(P_1))

write_csv(all_xlink, "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/ppi_scores_w-xlink_status_09182022.csv")
