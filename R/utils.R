col_OI <- function(df, pattern) {

  stopifnot("Not containing column of interest" =
              grepl(pattern, colnames(df)) %>% any())

  grep(pattern, colnames(df), value = TRUE)
}

clean <- function(df) {

  stopifnot("Input data does not contain correct column names" =
              c("Only.identified.by.site",
                "Reverse",
                "Potential.contaminant",
                "Gene.names",
                "Protein.IDs",
                "id") %in%
              colnames(df))

  df %>%
    dplyr::filter(.data$Only.identified.by.site != "+" &
                    .data$Reverse != "+" &
                    .data$Potential.contaminant != "+") %>%
    dplyr::rename(geneName = .data$Gene.names,
                  proteinID = .data$Protein.IDs) %>%
    tidyr::separate(.data$proteinID,
                    into = c("proteinID"),
                    sep = ";")
}

log2Signal <- function(df, col) {

  #log2 transformation on col
  df[ ,col] <- log2(df[ ,col])
  is.na(df[ ,col]) <- sapply(df[ ,col], is.infinite)

  df %>%
    dplyr::select(.data$id, .data$proteinID, .data$geneName, col)
}

longer <- function(df, col) {

  stopifnot("Cannot make longer: data does not contain correct intensity columns" =
              col %in% colnames(df))

  #aggregate experimental conditions into one col
  df %>%
    tidyr::pivot_longer(col,
                        names_to = "run",
                        values_to = "intensity") %>%
    tidyr::separate(.data$run,
                    into = c("experiment", "rep"),
                    sep = -1,
                    convert = TRUE)
}

wider <- function(df) {

  stopifnot("Cannot make wider: data does not contain correct columns" =
              c("experiment", "rep", "intensity") %in% colnames(df))

  df %>%
    tidyr::unite("run", .data$experiment, .data$rep,
                 sep = "") %>%
    tidyr::pivot_wider(names_from = .data$run,
                       values_from = .data$intensity)
}

selectValid <- function(df){

  df %>%
    dplyr::group_by(.data$id, .data$experiment) %>%
    dplyr::summarize(valid = sum(!is.na(.data$intensity))) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarize(maxValid = max(.data$valid)) %>%
    dplyr::filter(.data$maxValid >= 3) %>%
    dplyr::select(.data$id) %>%
    dplyr::left_join(df, by = "id")
}

impute2 <- function(df, LOD){

  stopifnot("Input data does not contain correct column names" =
              c("experiment",
                "intensity",
                "rep") %in%
              colnames(df))

  smy <- df %>%
    dplyr::group_by(.data$experiment) %>%
    dplyr::summarise(
      valid = sum(!is.na(.data$intensity)),
      var = stats::var(.data$intensity, na.rm = TRUE),
      mean = mean(.data$intensity, na.rm = TRUE)
    )

  stopifnot("At least one triplicates are all detected" =
              max(smy$valid) >= 3,
            "there is only one experimental group" =
              smy$experiment %>%
              unique() %>%
              length() > 1)

  sd_impute <- (smy %>%
                  dplyr::filter(.data$valid == max(smy$valid)))$var %>%
    min() %>% sqrt()

  set.seed(1)

  df %>%
    split(~experiment) %>%
    purrr::map(imputeMissing, LOD, sd_impute) %>%
    purrr::reduce(rbind)
}

imputeMissing <- function(df, mean_impute, sd_impute){

  valid <- (!is.na(df$intensity)) %>% sum()

  df[is.na(df$intensity),]$intensity <-
    switch(as.character(valid),
           "0" = stats::rnorm(is.na(df$intensity) %>% sum(),
                              mean = mean_impute,
                              sd = sd_impute),
           "1" = stats::rnorm(is.na(df$intensity) %>% sum(),
                              mean = mean(df$intensity, na.rm = TRUE),
                              sd = sd_impute),
           stats::rnorm(is.na(df$intensity) %>% sum(),
                        mean = mean(df$intensity, na.rm = TRUE),
                        sd = stats::var(df$intensity, na.rm = TRUE) ^ (1/2))
    )
  df
}

pFC <- function(df, h0, h1) {

  stopifnot("Input data does not contain correct column names" =
              c("experiment",
                "intensity",
                "rep",
                "id",
                "proteinID") %in%
              colnames(df))

  grp1 <- paste0("impute_LFQ.intensity.", h0)
  grp2 <- paste0("impute_LFQ.intensity.", h1)

  d1 <- df %>%
    dplyr::filter(.data$experiment == grp1)
  d2 <- df %>%
    dplyr::filter(.data$experiment == grp2)

  t <- stats::t.test(d2$intensity,
                     d1$intensity,
                     var.equal = TRUE)

  data.frame(id = unique(df$id),
             pvalue = t$p.value,
             log1P = -log10(t$p.value),
             log2fc = mean(d2$intensity) - mean(d1$intensity),
             t.value = t$statistic)
}
