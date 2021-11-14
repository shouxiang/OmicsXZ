library(tidyverse)
rm(list = ls())

mouse <- read.csv("D:/R program/proteomics/20210129 Tunicamycin MG132/proteinGroups/mouseUniprot.csv",
                  header = TRUE,
                  stringsAsFactors = FALSE)

mousePPI <- OmicsXZ::proteinDegree(mouse$Entry, species = 10090)

usethis::use_data(mousePPI, overwrite = TRUE)
