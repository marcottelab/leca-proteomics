library(tidyverse)

clstfile <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/results_lj/clustering.csv"
elutfile <- "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/results/leca_euks_elut.filtdollo.filt150p.raw.csv"
phylofile <- "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/ids_exps_species.raw.150p.p3.clade_filt.csv"

clst <- readr::read_csv(clstfile)
elut <- readr::read_csv(elutfile)
phy <- readr::read_csv(phylofile)

phy_fmt <- phy %>%
  select(ID, animals, plants, tsar, excavate) %>%
  separate(ID, into = c("ID1", "ID2"), sep = " ") %>%
  pivot_longer(ID1:ID2)

phy_fmt_uniq <- phy_fmt %>%
  select(-name) %>%
  rename(ID = value) %>% 
  unique()

phy_fmt_uniq[phy_fmt_uniq == 0] <- NA

phy_fmt_fill <- phy_fmt_uniq %>% 
  group_by(ID) %>%
  fill(everything(), .direction = "downup") %>%
  slice(1)

phy_fmt_fill[is.na(phy_fmt_fill)] <- 0
write_csv(phy_fmt_fill, "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/ids_clade_presence_absence.csv")

clst_fmt <- clst %>%
  select(ID, matches("cut*"), matches("human*"), matches("arath*")) %>%
  left_join(phy_fmt_fill, by = c("ID")) %>%
  select(ID, matches("cut_*"), animals, plants, tsar, excavate,
         matches("human*"), matches("arath*")) %>%
  left_join(elut, by = c("ID" = "X1"))

write_csv(clst_fmt, "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/results_lj/clustering_fmt.csv")