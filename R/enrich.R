#' Enrichment analysis including GO, KEGG, Reactome
#'
#' @param poi vector
#' @param u vector
#' @param Org string
#' @return list
#' @export
enrichAll <- function(poi, u, Org){

  e1 <- clusterProfiler::enrichGO(gene = poi,
                                  universe = u,
                                  OrgDb = switch(Org,
                                                 "mouse" = org.Mm.eg.db::org.Mm.eg.db,
                                                 "human" = org.Hs.eg.db::org.Hs.eg.db,
                                                 {
                                                   stop("Organism must be mouse or human")
                                                 }),
                                  keyType = "UNIPROT",
                                  ont = "CC",
                                  readable = TRUE) %>%
    clusterProfiler::simplify()

  CC <- e1@result %>%
    dplyr::mutate(Cat = "CC")

  e2 <- clusterProfiler::enrichGO(gene = poi,
                                  universe = u,
                                  OrgDb = switch(Org,
                                                 "mouse" = org.Mm.eg.db::org.Mm.eg.db,
                                                 "human" = org.Hs.eg.db::org.Hs.eg.db,
                                                 {
                                                   stop("Organism must be mouse or human")
                                                 }),
                                  keyType = "UNIPROT",
                                  ont = "BP",
                                  readable = TRUE) %>%
    clusterProfiler::simplify()
  BP <- e2@result %>%
    dplyr::mutate(Cat = "BP")

  e3 <- clusterProfiler::enrichGO(gene = poi,
                                  universe = u,
                                  OrgDb = switch(Org,
                                                 "mouse" = org.Mm.eg.db::org.Mm.eg.db,
                                                 "human" = org.Hs.eg.db::org.Hs.eg.db,
                                                 {
                                                   stop("Organism must be mouse or human")
                                                 }),
                                  keyType = "UNIPROT",
                                  ont = "MF",
                                  readable = TRUE) %>%
    clusterProfiler::simplify()
  MF <- e3@result %>%
    dplyr::mutate(Cat = "MF")

  #Convert UniProt Accession to ENTREZ ID
  x <- clusterProfiler::bitr(poi,
                             fromType = "UNIPROT", toType="ENTREZID",
                             OrgDb = switch(Org,
                                            "mouse" = "org.Mm.eg.db",
                                            "human" = "org.Hs.eg.db",
                                            {
                                              stop("Organism must be mouse or human")
                                            }))
  x <- x$ENTREZID
  b <- clusterProfiler::bitr(u,
                             fromType = "UNIPROT", toType="ENTREZID",
                             OrgDb = switch(Org,
                                            "mouse" = "org.Mm.eg.db",
                                            "human" = "org.Hs.eg.db",
                                            {
                                              stop("Organism must be mouse or human")
                                            }))
  b <- b$ENTREZID

  e4 <- clusterProfiler::enrichKEGG(gene = x,
                                    universe = b,
                                    organism = switch(Org,
                                                      "mouse" = "mmu",
                                                      "human" = "hsa",
                                                      {
                                                        stop("Organism must be mouse or human")
                                                      }))
  KEGG <- e4@result %>%
    dplyr::mutate(Cat = "KEGG")

  e5 <- ReactomePA::enrichPathway(gene = x,
                                  universe = b,
                                  organism = switch(Org,
                                                    "mouse" = "mouse",
                                                    "human" = "human",
                                                    {
                                                      stop("Organism must be mouse or human")
                                                    }))

  R <- e5@result %>%
    dplyr::mutate(Cat = "Reactome")

  result <- rbind(BP, CC) %>%
    rbind(MF) %>%
    rbind(KEGG) %>%
    rbind(R)

  result %>%
    tidyr::separate(.data$GeneRatio, into = c("poi_Cat", "poi_Qty"),
             sep = "/", convert = TRUE) %>%
    tidyr::separate(.data$BgRatio, into = c("bg_Cat", "bg_Qty"),
             sep = "/", convert = TRUE ) %>%
    dplyr::mutate(foldEnrichment = .data$poi_Cat/.data$poi_Qty/(.data$bg_Cat/.data$bg_Qty)) %>%
    dplyr::filter(.data$p.adjust < .05 & .data$poi_Cat >= 5) %>%
    dplyr::filter(.data$foldEnrichment > 3)
}
