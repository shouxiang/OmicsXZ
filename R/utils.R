#Names of all columns of signal intensity
getSignalColumnName <- function(df, pattern) {

  stopifnot("Not containing column of interest" =
              grepl(pattern, colnames(df)) %>% any())

  grep(pattern, colnames(df), value = TRUE)
}

clean <- function(df, col) {

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
    dplyr::rename(
      proteinID = .data$Protein.IDs) %>%
    tidyr::separate(.data$proteinID,
                    into = c("proteinID"),
                    sep = ";") %>%
    dplyr::select(.data$proteinID, col)
}

longer <- function(df, col, seperator) {

  stopifnot("Cannot make longer: data does not contain correct intensity columns" =
              col %in% colnames(df))

  if(seperator == "-1"){
    seperator <- -1
  }

  #aggregate experimental conditions into one col
  df %>%
    tidyr::pivot_longer(col,
                        names_to = "run",
                        values_to = "raw") %>%
    tidyr::separate(.data$run,
                    into = c("experiment", "rep"),
                    sep = seperator,
                    convert = TRUE)
}

log2Signal <- function(df) {

  stopifnot("Cannot make longer: data does not contain correct intensity columns" =
              "raw" %in% colnames(df))

  df[ ,"raw"] <- log2(df[ ,"raw"])

  df[is.infinite(df$raw),"raw"] <- NA

  df
}

selectValid <- function(df){

  df %>%
    dplyr::group_by(.data$proteinID, .data$experiment) %>%
    dplyr::summarize(valid = sum(!is.na(.data$raw))) %>%
    dplyr::group_by(.data$proteinID) %>%
    dplyr::summarize(maxValid = max(.data$valid)) %>%
    dplyr::filter(.data$maxValid >= 3) %>%
    dplyr::select(.data$proteinID)
}

impute <- function(df, width = 0.3,  downshift = 2.2) {

  stopifnot("Cannot make longer: data does not contain correct intensity columns" =
              "raw" %in% colnames(df))

  df <- df %>%
    dplyr::mutate(imputed  = raw)

  tmp.sd = width * stats::sd(df$raw, na.rm = TRUE)   # shrink sd width

  tmp.mean = mean(df$raw, na.rm = TRUE) -
    downshift * stats::sd(df$raw, na.rm = TRUE)   # shift mean of imputed values

  set.seed(1)

  df[is.na(df$imputed),"imputed"] =
    stats::rnorm(sum(is.na(df$imputed)), mean = tmp.mean, sd = tmp.sd)

  return(df)
}

get_P_foldChange <- function(df, h0, h1) {

  stopifnot("Input data does not contain correct column names" =
              c("proteinID",
                "experiment",
                "rep",
                "imputed") %in% colnames(df))

  grp1 <- paste0("LFQ.intensity.", h0)
  grp2 <- paste0("LFQ.intensity.", h1)

  d1 <- df %>%
    dplyr::filter(.data$experiment == grp1)
  d2 <- df %>%
    dplyr::filter(.data$experiment == grp2)

  t.result <- stats::t.test(d2$imputed,
                            d1$imputed,
                            var.equal = TRUE)

  data.frame(proteinID = unique(df$proteinID),
             pvalue = t.result$p.value,
             log1P = -log10(t.result$p.value),
             log2fc = mean(d2$imputed) - mean(d1$imputed),
             t.value = t.result$statistic,
             p.adjust = stats::p.adjust(t.result$p.value, method = "BH"))
}

wider <- function(df, values_from) {

  stopifnot("Cannot make wider: data does not contain correct columns" =
              c("proteinID","experiment", "rep", values_from) %in% colnames(df))

  df %>%
    dplyr::select(.data$proteinID, .data$experiment, .data$rep, !!dplyr::quo_name(values_from)) %>%
    tidyr::unite("run", .data$experiment, .data$rep,sep = "") %>%
    tidyr::pivot_wider(names_from = .data$run,
                       values_from = !!dplyr::quo_name(values_from))
}
