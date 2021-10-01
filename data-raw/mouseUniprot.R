library(tidyverse)

mouseUniprot <- read.csv("D:/R program/proteomics/20210129 Tunicamycin MG132/proteinGroups/mouseUniprot.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  mutate(Cys = str_count(Sequence, "C")) %>%
  rename(proteinID = Entry) %>%
  select(proteinID, Cys)

usethis::use_data(mouseUniprot, internal = TRUE)
