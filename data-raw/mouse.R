mouse <-
  mouseUniprot_site_location %>%
  left_join(mouseIUPred, by = "proteinID") %>%
  left_join(mousePPI, by = c("proteinID" = "UniProt"))

usethis::use_data(mouse, overwrite = TRUE)
