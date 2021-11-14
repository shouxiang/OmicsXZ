library(tidyverse)
rm(list = ls())

human <- read.csv("D:/R program/proteomics/20211018 PD/proteinGroups/humanUniprot.csv",
                  header = TRUE,
                  stringsAsFactors = FALSE)

humanPPI <- OmicsXZ::proteinDegree(human$Entry, species = 9606)

usethis::use_data(humanPPI, overwrite = TRUE)
