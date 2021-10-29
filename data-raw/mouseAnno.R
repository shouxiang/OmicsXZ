library(tidyverse)

mouse0 <- read.csv("D:/R program/proteomics/20210129 Tunicamycin MG132/proteinGroups/mouseUniprot.csv",
                   header = TRUE,
                   stringsAsFactors = FALSE)

mouseAnno <- mouse0 %>%
  OmicsXZ::celLoc() %>%
  OmicsXZ::cysFun() %>%
  select(proteinID, Cys, Nucleus, Membrane, Cytoplasm, Secreted,
         Mitochondrion, Golgi, ER, Endosome, Lysosome, Peroxisome,
         redoxActive, nucleophileCys, activeCys)

usethis::use_data(mouseAnno, overwrite = TRUE)
