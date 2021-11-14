library(tidyverse)

OmicsXZ::splitFasta(x = 13)

#select resultFasta dir
mouseUnfold <- OmicsXZ::mergeFasta()

usethis::use_data(mouseUnfold, overwrite = TRUE)
