library(tidyverse)

# cellular localization
human0 <- read.csv("D:/R program/proteomics/20211018 PD/proteinGroups/humanUniprot.csv",
                   header = TRUE,
                   stringsAsFactors = FALSE)

humanAnno <- human0 %>%
  OmicsXZ::celLoc() %>%
  OmicsXZ::cysFun() %>%
  select(proteinID, Cys, Nucleus, Membrane, Cytoplasm, Secreted,
         Mitochondrion, Golgi, ER, Endosome, Lysosome, Peroxisome,
         redoxActive, nucleophileCys, activeCys)

usethis::use_data(humanAnno, overwrite = TRUE)
