#' calculate the degrees of proteins
#'
#' @param proteinID char
#' @param species int
#'
#' @return dataframe
#' @export
proteinDegree <- function(proteinID, species) {

  stopifnot("species should not be char, no quotes" = is.double(species))

  string_bind <- STRINGdb::STRINGdb$new(version = "11",
                                        species = species,
                                        score_threshold = 400)

  string_id <- string_bind$map(data.frame(proteinID = proteinID),
                               "proteinID",
                               removeUnmappedRows = TRUE)

  string_id_agg <- string_id %>%
    dplyr::group_by(proteinID) %>%
    dplyr::summarise(
      stringID = list(.data$STRING_id)
    )

  degrees <- vector(length = nrow(string_id_agg))

  for (i in 1:nrow(string_id_agg)) {

    id <- string_id_agg[i, ]$stringID %>% unlist()
    x <- string_bind$get_neighbors(id)
    degrees[i] <- length(unique(intersect(x, string_id$STRING_id)))
  }

  string_id_agg %>%
    dplyr::select(-.data$stringID) %>%
    cbind(degrees)
}
