#' Volcano plots
#'
#' @param df data.frame
#' @param fcTH double
#' @return ggplot2.object
#' @importFrom ggplot2 ggplot geom_point geom_segment aes element_rect element_text margin theme theme_classic unit xlab ylab
#' @export

plotVolcano <- function(df, fcTH = 1.5){

  #data must contain required columns
  stopifnot("Input data does not contain correct column names" =
              c("log2fc",
                "pvalue",
                "p.adjust",
                "Cys",
                "log1P",
                "geneName") %in%
              colnames(df))

  xMin <- min(df$log2fc) - 1
  xMax <- max(df$log2fc) + 1
  yMax <- max(-log10(df$pvalue)) + 1
  pvalueTH <- max((df %>%
    dplyr::filter(.data$p.adjust < 0.05))$pvalue)

  hit <- df %>%
    dplyr::filter(.data$log2fc >= log2(fcTH) &
                    .data$pvalue <= pvalueTH &
                    .data$Cys > 0)

  ggplot(df, aes(.data$log2fc, .data$log1P)) +
    geom_point(colour = "grey", alpha = .4,
               size = 4) +
    geom_point(data = hit, colour = "red",
               size = 4) +
    geom_segment(aes(x = log2(fcTH), y = -log10(pvalueTH),
                     xend = xMax, yend = -log10(pvalueTH)),
                 linetype = "dashed", size = 1) +
    geom_segment(aes(x = log2(fcTH), y = -log10(pvalueTH),
                     xend = log2(fcTH), yend = yMax),
                 linetype = "dashed", size = 1) +
    theme_classic() +
    theme(
      strip.background = element_rect(colour = "white"),
      axis.text.x = element_text(margin = margin(t = 5),
                                 size = 40,
                                 face = "bold",
                                 color = "black"),
      axis.text.y = element_text(margin = margin(r = 5),size = 40,face = "bold",color = "black"),
      axis.ticks.length = unit(0.1, "in")) +
    xlab("") +
    ylab("")
  }
