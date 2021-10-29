#' Process proteomics data, two experimental conditions
#'
#' @param df data.frame
#' @param h0 string
#' @param h1 string
#' @param replicate int
#' @return data.frame
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data

processLFQ <- function(df, h0, h1, replicate = 3) {

  stopifnot("Iput experimental conditions are wrong" =
              paste0("LFQ.intensity.", rep(c(h0, h1), each = replicate), 1:replicate)
            %in% colnames(df) &
              h0 != h1)

  col <- col_OI(df, "^LFQ.intensity.")

  df1 <- preprocess(df, col)

  # imputation
  df2 <- df1 %>%
    longer(col)

  LOD <- min(df2$intensity, na.rm = TRUE) - 2

  df2 <- df2 %>%
    split(~id) %>%
    purrr::map(impute2, LOD) %>%
    purrr::reduce(rbind) %>%
    wider()

  colnames(df2) <- ifelse(colnames(df2) %in% col,
                          paste0("impute_", colnames(df2)),
                          colnames(df2))

  # calculate pvalue and fold change

  col <- col_OI(df2, "^impute_LFQ.intensity.")

  df3 <- df2 %>%
    longer(col) %>%
    split(~id) %>%
    purrr::map(pFC, h0, h1) %>%
    purrr::reduce(rbind) %>%
    dplyr::left_join(df1, by = "id") %>%
    dplyr::left_join(df2 %>%
                       dplyr::select(-.data$proteinID, -.data$geneName), by = "id") %>%
    dplyr::mutate(p.adjust = stats::p.adjust(.data$pvalue, method = "BH"))

  df3
}

preprocess <- function(df, col) {

  df %>%
    clean() %>%
    log2Signal(col) %>%
    longer(col) %>%
    #select intensity detected in all triplicates in
    #at least one experimental group
    selectValid() %>%
    wider()
}
