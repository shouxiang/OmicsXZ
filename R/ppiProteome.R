
#' export the degree and betweenness of the proteome
#'
#' @param species string
#' @param score_threshold int
#'
#' @return data.frame
#' @export
#'
ppiProteome <- function(species, score_threshold = 600) {

  if(species == "mouse"){

    organismID <- 10090

    proteome <-
      utils::read.csv(
        "D:/R program/proteomics/database/mouseUniprot_site_location.csv",
        header = TRUE,
        stringsAsFactors = FALSE) %>%
      dplyr::select(.data$Entry)
  } else if(species == "human") {

    organismID <- 9606

    proteome <-
      utils::read.csv(
        "D:/R program/proteomics/database/humanUniprot_site_location.csv",
        header = TRUE,
        stringsAsFactors = FALSE) %>%
      dplyr::select(.data$Entry)
  } else {
    stop("only support mouse and human.")
  }

  string_db <-
    STRINGdb::STRINGdb$new(version = "11",
                           species = organismID,
                           score_threshold = score_threshold)

  uniprot_Stringid <-
    string_db$map(
      data.frame(UniProt = proteome$Entry),
      "UniProt",
      removeUnmappedRows = TRUE)

  g0 <- string_db$get_graph()

  g0_edgeList <-
    igraph::as_edgelist(g0) %>%
    unique()

  g <-
    igraph::graph_from_edgelist(
      g0_edgeList,
      directed = FALSE
    )

  g_degree <- igraph::degree(g)

  g_between <-
    igraph::betweenness(
      g,
      directed = FALSE)

  degree.df <-
    tibble::tibble(STRING_id = names(g_degree), degree0 = g_degree)
  betweenness.df <-
    tibble::tibble(STRING_id = names(g_between), betweenness0 = g_between)

  result <-
    dplyr::inner_join(
      uniprot_Stringid, degree.df, by = "STRING_id") %>%
    dplyr::left_join(
      betweenness.df, by = "STRING_id"
    ) %>%
    dplyr::group_by(.data$UniProt) %>%
    dplyr::summarise(
      degree = max(.data$degree0),
      betweenness = max(.data$betweenness0)
    )

  result
}
