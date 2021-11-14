library(tidyverse)

OmicsXZ::splitFasta(x = 15)

#select resultFasta dir
humanUnfold <- OmicsXZ::mergeFasta()

usethis::use_data(humanUnfold, overwrite = TRUE)
