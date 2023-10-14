human <-
  humanUniprot_site_location %>%
  left_join(humanIUPred, by = "proteinID") %>%
  left_join(humanPPI, by = c("proteinID" = "UniProt"))

usethis::use_data(human, overwrite = TRUE)
