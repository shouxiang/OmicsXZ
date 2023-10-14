#' 2 experimental conditions
#'
#' @param df data.frame
#' @param h0 string
#' @param h1 string
#' @param seperator string
#'
#' @return data.frame
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data

process2groups <- function(df, h0, h1, seperator = "-1") {

  df1 <- processLFQ(df, seperator)

  df_P_fc <- df1 %>%
    split(~proteinID) %>%
    purrr::map(get_P_foldChange, h0, h1) %>%
    purrr::reduce(rbind)

  #new column name for imputed value
  df2 <- df1 %>%
    wider("imputed")
  columnNames <- getSignalColumnName(df2, "^LFQ.intensity.")
  colnames(df2) <-
    ifelse(colnames(df2) %in% columnNames,
           paste0("impute_", colnames(df2)),
           colnames(df1))

  df3 <- df1 %>%
    wider("raw")

  df_all <-
    dplyr::left_join(df3,df2,by = "proteinID") %>%
    dplyr::left_join(df_P_fc, by = "proteinID")

  df_all
}

processLFQ <- function(df, seperator) {

  columnNames <- getSignalColumnName(df, "^LFQ.intensity.")

  df_log2 <- df %>%
    clean(columnNames) %>%
    longer(columnNames, seperator) %>%
    log2Signal()

  df_valid <- df_log2 %>%
    selectValid() %>%
    dplyr::left_join(df_log2, by = "proteinID")

  df_imputed <- df_valid %>%
    split(~experiment) %>%
    purrr::map(impute) %>%
    purrr::reduce(rbind)

  df_imputed
}

