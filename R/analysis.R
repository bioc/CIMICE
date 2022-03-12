# This script contains functions to quicken the
# data analysis and to prepare working examples
utils::globalVariables(c("A","B","C","D","freq"))

#' Creates a simple example dataset
#'
#' @return a simple mutational matrix
#'
#' @examples
#' example_dataset()
#'
#' @export example_dataset
example_dataset <- function(){
    make_dataset(A,B,C,D) %>% # genes
        # samples
        update_df("S1", 0, 0, 0, 1) %>%
        update_df("S2", 1, 0, 0, 0) %>%
        update_df("S3", 1, 0, 0, 0) %>%
        update_df("S4", 1, 0, 0, 1) %>%
        update_df("S5", 1, 1, 0, 1) %>%
        update_df("S6", 1, 1, 0, 1) %>%
        update_df("S7", 1, 0, 1, 1) %>%
        update_df("S8", 1, 1, 0, 1)
}

#' Creates a simple example dataset with frequency column
#'
#' @return a simple mutational matrix
#'
#' @examples
#' example_dataset_withFreqs()
#'
#' @export example_dataset_withFreqs
example_dataset_withFreqs <- function(){
    dataset <- make_dataset(A,B,C,D) %>% # genes
        # samples
        update_df("G1", 0, 0, 0, 1) %>%
        update_df("G2", 1, 0, 0, 0) %>%
        update_df("G3", 1, 0, 0, 1) %>%
        update_df("G4", 1, 1, 0, 1) %>%
        update_df("G5", 1, 0, 1, 1) %>%
        update_df("G6", 1, 1, 0, 1)
    counts <- c(1,2,1,2,1,1)
    list(matrix = dataset, counts = counts)
}

#' Run CIMICE preprocessing
#'
#' executes the preprocessing steps of CIMICE
#'
#' Preprocessing steps:
#'
#' 1) dataset is compacted
#'
#' 2) genotype frequencies are computed
#'
#' 3) labels are prepared
#'
#' @param dataset a mutational matrix as a (sparse) matrix
#'
#' @return a list containing the mutational matrix ("samples"),
#' the mutational frequencies of the genotypes ("freqs"),
#' the node labels ("labels") and finally the gene names ("genes")
#'
#' @examples
#' require(dplyr)
#' example_dataset() %>% dataset_preprocessing
#'
#' @export dataset_preprocessing
dataset_preprocessing <- function(dataset){
    # compact
    compactedDataset <- compact_dataset(dataset)
    # run preprocessing on genotypes and frequencies
    dataset_preprocessing_population(compactedDataset)
}

#' Run CIMICE preprocessing for poulation format dataset
#'
#' executes the preprocessing steps of CIMICE
#'
#' Preprocessing steps:
#'
#' 1) genotype frequencies are computed
#'
#' 2) labels are prepared
#'
#' @param compactedDataset a list (matrix: a mutational matrix, 
#' counts: number of samples with given genotype). 
#' "counts" is normalized automatically. 
#'
#' @return a list containing the mutational matrix ("samples"),
#' the mutational frequencies of the genotypes ("freqs"),
#' the node labels ("labels") and finally the gene names ("genes")
#'
#' @examples
#' require(dplyr)
#' example_dataset_withFreqs() %>% dataset_preprocessing_population
#'
#' @export dataset_preprocessing_population
dataset_preprocessing_population <- function(compactedDataset){
    # dataset is already compacted by hypothesys
    # prepare for the analysis
    samples <- compactedDataset$matrix
    # save genes' names
    genes <- colnames(compactedDataset$matrix)
    # keep the information on frequencies for further analysis
    freqs <- compactedDataset$counts/sum(compactedDataset$counts)
    # prepare node labels listing the mutated genes for each node
    labels <- prepare_labels(samples, genes)
    # prepare node labels displaying genotype names
    if( is.null(compactedDataset$row_names) ){
        compactedDataset$row_names <- rownames(compactedDataset$matrix)
    }
    matching_samples <- compactedDataset$row_names
    # fix Colonal genotype absence, if needed
    fix <- fix_clonal_genotype(samples, freqs, labels, matching_samples)
    samples <- fix[["samples"]]
    freqs <- fix[["freqs"]]
    labels <- fix[["labels"]]
    matching_samples <- fix[["matching_samples"]]
    # return a list with the prepared dataset and its additional information
    list("samples" = samples, "freqs" = freqs,
        "labels" = labels, "genes" = genes, "matching_samples" = matching_samples)
}

#' Default preparation of graph topology
#'
#' By default, CIMICE computes the relation between
#' genotypes using the subset relation.
#' For the following steps it is also important
#' that the transitive edges are removed.
#'
#' @param samples mutational matrix
#' @param labels genotype labels
#'
#' @return a graph with the wanted topology
#'
#' @examples
#' require(dplyr)
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' graph_non_transitive_subset_topology(samples, labels)
#' @export graph_non_transitive_subset_topology
graph_non_transitive_subset_topology <- function(samples, labels){
    # compute edges based on subset relation
    edges <- build_topology_subset(samples)
    # remove transitive edges and prepare igraph object
    build_subset_graph(edges, labels)
}

#' Compute default weights
#'
#' This procedure computes the weights for edges of a
#' graph accordingly to CIMICE specification.
#' (See vignettes for further explainations)
#'
#' @param g a graph (must be a DAG with no transitive edges)
#' @param freqs observed frequencies of genotypes
#'
#' @return a graph with the computed weights
#'
#' @examples
#' require(dplyr)
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' compute_weights_default(g, freqs)
#'
#' @export compute_weights_default
compute_weights_default <- function(g, freqs){
    # prepare adj matrix
    A <- as_adj(g)
    # pre-compute exiting edges from each node
    no.of.children <- get_no_of_children(A,g)
    # compute the four steps
    upWeights <- computeUPW(g, freqs, no.of.children, A)
    normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
    downWeights <- computeDWNW(g, freqs, no.of.children, A, normUpWeights)
    normDownWeights <- normalizeDWNW(g, freqs, no.of.children, A, downWeights)
    normDownWeights
}

#' Run CIMICE defaults
#'
#' This function executes CIMICE analysis on a dataset using default settings.
#'
#' @param dataset a mutational matrix as a (sparse) matrix
#' @param mode indicates the used input format. Must be either "CAPRI" or "CAPRIpop"
#'
#' @return a list object representing the graph computed by CIMICE with the
#' structure `list(topology = g, weights = W, labels = labels, freqs=freqs)`
#'
#' @examples
#' quick_run(example_dataset())
#'
#' @export quick_run
quick_run <- function(dataset, mode="CAPRI"){
    # preprocess
    preproc <- NULL
    if(mode == "CAPRI"){
        preproc <- dataset_preprocessing(dataset)
    }else if(mode == "CAPRIpop"){
        preproc <- dataset_preprocessing_population(dataset)
    }else{
        stop("Unsupported input mode", mode, "use CAPRI o CAPRIpop")
    }
    samples <- preproc[["samples"]]
    freqs   <- preproc[["freqs"]]
    labels  <- preproc[["labels"]]
    genes   <- preproc[["genes"]]
    matching_samples   <- preproc[["matching_samples"]]
    g <- graph_non_transitive_subset_topology(samples,labels)
    W <- compute_weights_default(g, freqs)
    list(topology = g, weights = W, labels = labels, freqs=freqs, matching_samples=matching_samples)
}



