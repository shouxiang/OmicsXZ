#' calculate percentages, fold enrichment and significance of an annotation
#'
#' @param anno string
#' @param df data.frame
#' @param ref data.frame
#'
#' @return data.frame
#' @export

annoInfo <- function(anno, df, ref) {

  stopifnot("dataframe of interest does not contain functional annotations" =
              c(anno, "proteinID") %in%
              colnames(df),
            "reference does not contain functional annotations" =
              c(anno, "proteinID") %in%
              colnames(ref))

  poi_Cat <- sum(df[[anno]])
  poi_NotCat <- nrow(df) - poi_Cat

  NotPoi <- ref %>%
    dplyr::filter(!.data$proteinID %in% df$proteinID)
  NotPoi_Cat <- sum(NotPoi[[anno]])
  NotPoi_NotCat <- nrow(NotPoi) - NotPoi_Cat

  p.value <- (data.frame(c(poi_Cat, poi_NotCat),
                         c(NotPoi_Cat, NotPoi_NotCat)) %>%
                stats::fisher.test())$p.value

  sig <- cut(p.value,
             breaks = c(0, .0001, .001, .01, .05, Inf),
             labels = c("****",
                        "***",
                        "**",
                        "*",
                        "ns"),
             right = FALSE) %>% as.character()

  foldEnrichment = sum(df[[anno]])/nrow(df)/
    (sum(ref[[anno]])/nrow(ref))

  data.frame(annotation = anno,
             data = deparse(substitute(df)),
             percentage = mean(df[[anno]]),
             foldEnrichment = foldEnrichment,
             p.value = p.value,
             label = paste0(round(foldEnrichment, digits = 2), ", ", sig))
}

celLoc <- function(df) {

  stopifnot("Input data does not contain correct column names" =
              c("Sequence",
                "Subcellular.location..CC.",
                "Entry"
                ) %in%
              colnames(df))

  df %>%
    dplyr::mutate(Cys = stringr::str_count(.data$Sequence, "C"),
                  Nucleus = grepl("Nucleus|Chromosome", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Membrane = grepl("Membrane", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Cytoplasm = grepl("Cytoplasm", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Secreted = grepl("Secreted", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Mitochondrion = grepl("Mitochondrion", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Golgi = grepl("Golgi apparatus", .data$Subcellular.location..CC., ignore.case = TRUE),
                  ER = grepl("Endoplasmic reticulum", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Endosome = grepl("Endosome", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Lysosome = grepl("Lysosome", .data$Subcellular.location..CC., ignore.case = TRUE),
                  Peroxisome = grepl("Peroxisome", .data$Subcellular.location..CC., ignore.case = TRUE)) %>%
    dplyr::rename(proteinID = .data$Entry)
}

# Cys functional annotations
cysFun <- function(df){

  stopifnot("Input data does not contain Cys functional annotations" =
              c("Sequence",
                "Active.site",
                "Disulfide.bond"
              ) %in%
              colnames(df))

  nucleophileCys <- rep(FALSE, nrow(df))
  activeCys <- rep(FALSE, nrow(df))

  for(i in 1:nrow(df)){

    nucleophileCys[i] <- isNucleophileCys(df[i,]$Active.site, df[i,]$Sequence)
    activeCys[i] <- isActiveCys(df[i,]$Active.site, df[i,]$Sequence)
  }

  df %>%
    cbind(nucleophileCys) %>%
    cbind(activeCys) %>%
    dplyr::mutate(redoxActive = stringr::str_detect(.data$Disulfide.bond, "Redox-active"))
}

isNucleophileCys <- function(annotation, aa){

  result <- FALSE
  str <- stringr::str_extract_all(annotation,
                                  "ACT_SITE \\d*;  /note=\"Nucleophile\"")
  str[lengths(str) == 0] <- NA

  for(i in 1:length(str[[1]])){

    if(!is.na(str[[1]][i])){
      cysSite <- stringr::str_remove_all(str[[1]][i],
                                         "(ACT_SITE |;  /note=\"Nucleophile\")") %>%
        as.integer()
      if(substr(aa, cysSite, cysSite) == "C"){
        result <- TRUE
        break
      }
    }
  }
  result
}

isActiveCys <- function(annotation, aa){

  result <- FALSE
  str <- stringr::str_extract_all(annotation,
                                  "ACT_SITE \\d*;")
  str[lengths(str) == 0] <- NA

  for(i in 1:length(str[[1]])){

    if(!is.na(str[[1]][i])){
      cysSite <- stringr::str_remove_all(str[[1]][i],
                                         "ACT_SITE |;") %>%
        as.integer()
      if(substr(aa, cysSite, cysSite) == "C"){
        result <- TRUE
        break
      }
    }
  }
  result
}
