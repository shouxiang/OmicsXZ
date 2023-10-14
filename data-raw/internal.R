library(tidyverse)
rm(list = ls())

# mouse -------------------------------------------------------------------

mouseUniprot_site_location.df <-
  read.csv(
    "D:/Publication/TME/sourceData/ms/database/uniprot_mouse_cys_loc.csv",
    header = TRUE,
    stringsAsFactors = FALSE) %>%
  mutate(geneName = map_chr(
    Gene.Names,
    function(text){
      strsplit(text, " ")[[1]][1]
    }))

{## Functional Annotation --------------------------------------------------

  mouseUniprot_site_location <-
    mouseUniprot_site_location.df %>%
    OmicsXZ::celLoc() %>%
    OmicsXZ::cysFun() %>%
    select(Nucleus, Membrane, Cytoplasm, Secreted,
           Mitochondrion, Golgi, ER, Endosome, Lysosome, Peroxisome,
           redoxActive, nucleophileCys, activeCys,
           proteinID, geneName, Cys)
}

{## IUPred ------------------------------------------------------------------

  # OmicsXZ::splitFasta(x = 13)

  #select resultFasta dir
  mouseIUPred <-
    OmicsXZ::mergeFasta(
      "D:/R program/proteomics/20210129 Tunicamycin MG132/IUPred2A/resultFasta")
}

{## PPI ---------------------------------------------------------------------

  mousePPI <-
    OmicsXZ::ppiProteome("mouse")

}

# human -------------------------------------------------------------------

humanUniprot_site_location.df <-
  read.csv(
    "D:/Publication/TME/sourceData/ms/database/uniprot_human_cys_loc.csv",
    header = TRUE,
    stringsAsFactors = FALSE) %>%
  mutate(geneName = map_chr(
    Gene.Names,
    function(text){
      strsplit(text, " ")[[1]][1]
    }))

{## Functional Annotation --------------------------------------------------

  humanUniprot_site_location <-
    humanUniprot_site_location.df %>%
    OmicsXZ::celLoc() %>%
    OmicsXZ::cysFun() %>%
    select(Nucleus, Membrane, Cytoplasm, Secreted,
           Mitochondrion, Golgi, ER, Endosome, Lysosome, Peroxisome,
           redoxActive, nucleophileCys, activeCys,
           proteinID, geneName, Cys)
}

{## IUPred ------------------------------------------------------------------

  # OmicsXZ::splitFasta(x = 15)

  #select resultFasta dir
  humanIUPred <-
    OmicsXZ::mergeFasta(
      "D:/R program/proteomics/20211018 PD/IUPred2A/resultFasta")
}

{## PPI ---------------------------------------------------------------------

  humanPPI <-
    OmicsXZ::ppiProteome("human")

}

# Save --------------------------------------------------------------------

usethis::use_data(
  mouseUniprot_site_location,
  mouseIUPred,
  mousePPI,
  humanUniprot_site_location,
  humanIUPred,
  humanPPI,
  internal = TRUE,
  overwrite = TRUE)
