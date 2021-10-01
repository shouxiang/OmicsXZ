#' Preprocess proteomics data
#'
#' @param df dataframe
#' @param left string
#' @param right string
#'
#' @return dataframe
#' @importFrom magrittr %>%
#' @importFrom rlang .data
preprocess <- function(df, left, right) {
  #data must contain required columns
  stopifnot("Input data does not contain correct column names" =
              c("Only.identified.by.site",
                "Reverse",
                "Potential.contaminant",
                "Gene.names",
                "Protein.IDs") %in%
              colnames(df),
            "Iput experimental conditions are wrong" =
              paste0("LFQ.intensity.", rep(c(left, right), each = 3), 1:3)
            %in% colnames(df) &
              left != right)

  df <- df %>%
    dplyr::filter(.data$Only.identified.by.site != "+" &
                    .data$Reverse != "+" &
                    .data$Potential.contaminant != "+") %>%
    dplyr::rename(geneName = .data$Gene.names,
                  proteinID = .data$Protein.IDs) %>%
    tidyr::separate(.data$proteinID,
                    into = c("proteinID"),
                    sep = ";")

  tmp <- df[ ,grepl("^LFQ.intensity.", colnames(df))]

  colnames(tmp) <- stringr::str_remove(colnames(tmp), "^LFQ.intensity.")

  tmp <- log2(tmp)

  #log2 transformation
  df1 <- df %>%
    dplyr::select(.data$id, .data$proteinID, .data$geneName) %>%
    cbind(tmp)

  is.na(df1) <- sapply(df1, is.infinite)

  #select intensity detected in all triplicates in
  #at least one experimental group
  df2 <- df1[is.na(df1[[paste0(left, 1)]]) +
               is.na(df1[[paste0(left, 2)]]) +
               is.na(df1[[paste0(left, 3)]]) == 0 |
               is.na(df1[[paste0(right, 1)]]) +
               is.na(df1[[paste0(right, 2)]]) +
               is.na(df1[[paste0(right, 3)]]) == 0,]

  df3 <- df2 %>%
    tidyr::pivot_longer(c(paste0(left, 1:3), paste0(right, 1:3)),
                        names_to = "run",
                        values_to = "intensity") %>%
    tidyr::separate(.data$run,
                    into = c("experiment", "rep"),
                    sep = nchar(left),
                    convert = TRUE)

  LOD <- min(df3$intensity, na.rm = TRUE) - 2

  df4 <- df3 %>%
    split(~id) %>%
    purrr::map(impute, LOD) %>%
    purrr::reduce(rbind)

  df5 <- df4 %>%
    tidyr::unite("run", .data$experiment, .data$rep,
                 sep = "")

  df5$run <- paste0(df5$run, "_impute")

  df6 <- df5 %>%
    tidyr::pivot_wider(names_from = .data$run,
                       values_from = .data$intensity)

  df7 <- df4 %>%
    split(~id) %>%
    purrr::map(pFC, left, right) %>%
    purrr::reduce(rbind) %>%
    dplyr::left_join(df2, by = "id") %>%
    dplyr::left_join(df6 %>%
                       dplyr::select(-.data$proteinID, -.data$geneName), by = "id") %>%
    dplyr::left_join(mouseUniprot %>%
                       dplyr::select(.data$proteinID, .data$Cys), by = "proteinID")

  df7
}

#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats rnorm
#' @importFrom stats var
impute <- function(df, LOD) {

  stopifnot("Input data does not contain correct column names" =
              c("experiment",
                "intensity",
                "rep") %in%
              colnames(df))

  smy <- df %>%
    dplyr::group_by(.data$experiment) %>%
    dplyr::summarise(
      na = sum(is.na(.data$intensity)),
      s2 = stats::var(.data$intensity),
      mean = mean(.data$intensity, na.rm = TRUE)
    )

  #At least one group contains no missing values
  stopifnot("At least one triplicates are all detected" =
              min(smy$na) == 0)

  na <- max(smy$na)

  full <- smy %>%
    dplyr::filter(.data$na == 0)

  miss <- smy %>%
    dplyr::filter(.data$na != 0)

  set.seed(1)

  if(na == 3) {
    df[is.na(df$intensity),]$intensity <-
      stats::rnorm(miss$na,
                   mean = LOD,
                   sd = full$s2 ^ (1/2))
  }else {
    df[is.na(df$intensity),]$intensity <-
      stats::rnorm(miss$na,
                   mean = miss$mean,
                   sd = full$s2 ^ (1/2))
  }
  df
}

#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats t.test
pFC <- function(df, left, right) {

  stopifnot("Input data does not contain correct column names" =
              c("experiment",
                "intensity",
                "rep",
                "id",
                "proteinID") %in%
              colnames(df))

  d1 <- df %>%
    dplyr::filter(.data$experiment == left)
  d2 <- df %>%
    dplyr::filter(.data$experiment == right)

  t <- stats::t.test(d1$intensity,
                     d2$intensity,
                     var.equal = TRUE)

  data.frame(id = unique(df$id),
             pvalue = t$p.value,
             log1P = -log10(t$p.value),
             log2fc = mean(d1$intensity) - mean(d2$intensity),
             t.value = t$statistic)
}
