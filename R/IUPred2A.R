#' split fasta files into smaller ones
#'
#' @param inputFasta path
#' @param x int
#'
#' @return NULL
#' @export
splitFasta <- function(inputFasta = file.choose(), x) {

  stopifnot("input should be a fasta file" = grepl("\\.fasta$", inputFasta))

  fasta <- data.frame(txt = readLines(inputFasta)) %>%
    dplyr::mutate(grp = cumsum(grepl('^>sp\\|', .data$txt)))

  protein_count <- max(fasta$grp)

  t <- stringr::str_split(inputFasta, "\\\\", simplify = TRUE)
  dir <- paste0(t[1:length(t) - 1], collapse = "\\") %>%
    paste0("\\splitFasta")

  if(!file.exists(dir)){
    dir.create(dir)
  }

  j <- 0
  for(i in 1 : x){

    tmp <- fasta %>%
      dplyr::filter(.data$grp > as.integer((i-1)/x * protein_count) &
                      .data$grp <= as.integer(i/x * protein_count)
      )

    path <- paste0(dir, "\\", i, ".fasta")
    writeLines(tmp$txt, path)

    j <- j + max(tmp$grp) - min(tmp$grp) + 1
  }
}


#' Merge and summarize IUPred results
#'
#' @param inputDir path
#'
#' @return NULL
#' @export
mergeFasta <- function(inputDir) {

  file_name <- list.files(path = inputDir, pattern = '^\\d+.result$')

  file_list <- paste0(inputDir, "\\", file_name)

  purrr::map(file_list, readIUPred) %>%
    purrr::reduce(rbind)
}

readIUPred <- function(file) {

  data.frame(txt = readLines(file)) %>%
    dplyr::mutate(grp = cumsum(grepl('# IUPred2A.+', .data$txt))) %>%
    dplyr::group_by(.data$grp) %>%
    dplyr::mutate(proteinID = gsub(pattern = '^.*\\|(.{6,10})\\|.*',
                                   replacement = '\\1',
                                   x = grep('>sp', .data$txt, value = TRUE))) %>%
    dplyr::slice(6:(dplyr::n()-3)) %>%
    tidyr::separate(.data$txt, c("pos", "aa", "IUPred", "ANCHOR"), '\t', convert = TRUE) %>%
    dplyr::group_by(.data$proteinID) %>%
    dplyr::summarise(meanIUPred = mean(.data$IUPred),
                     meanANCHOR = mean(.data$ANCHOR))
}

