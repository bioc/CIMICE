#' CIMICE Package
#'
#' CIMICE-R: (Markov) Chain Method to Infer Cancer Evolution
#'
#' @description     R implementation of the CIMICE tool.
#' CIMICE is a tool in the field of tumor phylogenetics and
#' its goal is to build a Markov Chain (called Cancer Progression Markov Chain, CPMC) in order to model tumor subtypes evolution.
#' The input of CIMICE is a Mutational Matrix, so a boolean matrix representing altered genes in a
#' collection of samples. These samples are assumed to be obtained with single-cell DNA analysis techniques and
#' the tool is specifically written to use the peculiarities of this data for the CMPC construction.
#' See `https://github.com/redsnic/tumorEvolutionWithMarkovChains/tree/master/GenotypeEvolutionPaths` for the
#' original Java version of this tool.
#'
#' @author Nicolò Rossi \email{olocin.issor@gmail.com}
#'
#' @docType package
#' @name CIMICE
#'
#' @importFrom dplyr select enexprs select_if %>% ungroup summarise group_by_all n mutate_all count group_by 
#' @importFrom dplyr inner_join left_join mutate pull
#' @importFrom ggplot2 ggplot geom_histogram aes arrow unit labs
#' @importFrom glue glue
#' @importFrom tidyr drop_na
#' @importFrom igraph V E as_ids tail_of head_of get.edge.ids
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
#' @importFrom igraph graph_from_edgelist as_adj E<- V<-
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom visNetwork toVisNetworkData visNetwork visHierarchicalLayout
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom purrr map2 map_dbl imap_chr map map_chr map_int map2_chr rbernoulli
#' @importFrom ggraph ggraph geom_node_point geom_node_label geom_edge_loop
#' @importFrom ggraph geom_edge_link label_rect
#' @importFrom stats cor
#' @importFrom utils read.csv
#' @importFrom maftools read.maf mutCountMatrix
#' @importFrom assertthat are_equal assert_that
#' @importFrom tidygraph %E>% %N>% activate as_tbl_graph bind_edges
#' @importFrom expm %^%
#' @import Matrix
#'
NULL
