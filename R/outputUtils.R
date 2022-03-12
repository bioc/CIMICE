utils::globalVariables(c("node1.label", "node2.label"))
#' ggplot graph output
#'
#' Draws the output graph using ggplot
#'
#' @param out the output object of CIMICE (es, from quick run)
#' @param digits precision for edges' weights
#' @param ... other arguments for format_labels
#'
#' @return ggraph object representing g as described
#'
#' @examples
#' draw_ggraph(quick_run(example_dataset()))
#' 
#' @export draw_ggraph
draw_ggraph <- function(out, digits = 4, ...){
    g <- out$topology 
    W <- out$weights
    labels <- make_labels(out, ...)
    
    V(g)$label <- labels

    weights <- signif(W, digits = 4)

    ggraph(g, layout = "sugiyama") +
        geom_node_point() +
        geom_node_label(aes(label=labels)) +
        geom_edge_link(
            aes(label=weights, start_cap = label_rect(node1.label),
                end_cap = label_rect(node2.label)),
            arrow = arrow(length = unit(4, 'mm')),
            label_dodge = unit(2.5, 'mm'),
            angle_calc = 'along'
        )
}

#' Format labels for output object
#'
#' Prepare labels based on multiple identifiers so that they do not excede
#' a certain size (if they do, a simple number is used)
#' 
#' @param labels a charachter vector of the labels to manage
#' @param max_col maximum number of identifiers in a single row for a label
#' @param max_row maximum number of rows of identifiers in a label
#' 
#' @return the updated labels
#' 
#' @examples  
#' format_labels(c("A, B", "C, D, E"))
#'
#' @export format_labels
format_labels <- function(labels, max_col = 3, max_row = 3){
    # prepare basic labels as replacement of too long labels
    basic_labels <- map_chr(seq(1,length(labels)), ~ as.character(.))
    for(i in seq(1,length(labels))){
        # divide the label in tokens
        lsplit <- strsplit(labels[i], ", ")[[1]]
        if(length(lsplit) > max_col*max_row){
            # label too long
            labels[i] <- basic_labels[i]
        }else{
            # creation of the label's grid of IDs
            new_label <- ""
            col <- 1
            first_run <- TRUE
            for(l in lsplit){
                
                if(col == 1){
                    if(col == max_col){
                        col <- 1
                    }else{
                        col <- col + 1    
                    }
                    if(first_run){
                        new_label <- c(new_label, l)
                        first_run <- FALSE
                    }else{
                        new_label <- c(new_label, "\n", l)
                    }
                }else if(col == max_col){
                    col <- 1
                    new_label <- c(new_label, ", ", l)
                }else{
                    col <- col + 1
                    new_label <- c(new_label,", ", l)
                }
            }
            labels[i] <- paste(new_label, collapse = "")
        }
    }
    labels
}

#' Helper function to create labels 
#' 
#' This function helps creating labels for nodes with different information
#' 
#' @param out the output object of CIMICE (es, from quick run)
#' @param mode which labels to print: samplesIDs (matching samples), sequentialIDs (just a sequential numer), geneIDs (genotype)
#' @param max_col identifiers are represented in a max_col times max_row grid (if the number of IDs exceeds, a the sequentialID number is used instead) 
#' @param max_row identifiers are represented in a max_col times max_row grid (if the number of IDs exceeds, a the sequentialID number is used instead)
#' 
#' @return the requested labels
#'
#' @examples  
#' make_labels(quick_run(example_dataset()))
#' 
#' @export make_labels
make_labels <- function(out, mode = "samplesIDs", max_col = 3, max_row = 3){
    if(mode == "samplesIDs"){
        labels <- out$matching_samples %>% as.character %>% format_labels(max_col, max_row)
    }else if(mode == "sequentialIDs"){
        labels <- map_chr(seq(1,length(out$labels)), ~ as.character(.)) 
    }else if(mode == "geneIDs"){
        labels <- out$labels %>% as.character %>% format_labels(max_col, max_row)
    }else{
        stop("Invalid output mode:", mode)
    }
    labels
}


#' NetworkD3 graph output
#'
#' Draws the output graph using networkD3
#'
#' @param out the output object of CIMICE (es, from quick run)
#' @param ... other arguments for format_labels
#'
#' @return networkD3 object representing g as described
#'
#' @examples
#' draw_networkD3(quick_run(example_dataset()))
#'
#' @export draw_networkD3
draw_networkD3 <- function(out, ...){
    g <- out$topology 
    W <- out$weights
    labels <- make_labels(out, ...)
    
    V(g)$label <- labels
    E(g)$label <- W

    gr <- igraph_to_networkD3(g)
    gr$nodes$label <- unlist(labels, use.names=FALSE)
    # Create force directed network plot
    p<-forceNetwork(Links = gr$links, Nodes = gr$nodes,
                    Source = 'source', Target = 'target', Value = "value",
                    NodeID = 'label', Group = 'name',
                    zoom=TRUE, arrows = TRUE, opacityNoHover = FALSE)
    p
}

# #' === DEPRECATED ===
# 
# VisNetwork graph output (default)
#
# Draws the output graph using VisNetwork
#
# @param g graph to be drawn
# @param W weights on edges
# @param labels node labels
# @param original_dataset if not NULL, recover matching samples from the rownames of the original dataset
#
# @return visNetwork object representing g as described
#
# @examples
# require(dplyr)
# preproc <- example_dataset() %>% dataset_preprocessing
# samples <- preproc[["samples"]]
# freqs   <- preproc[["freqs"]]
# labels  <- preproc[["labels"]]
# genes   <- preproc[["genes"]]
# g <- graph_non_transitive_subset_topology(samples, labels)
# W <- compute_weights_default(g, freqs)
# draw_visNetwork(g,W,labels)
#
# @export draw_visNetwork
# draw_visNetwork <- function(out, original_dataset = NULL){
# 
#     V(g)$label <- labels
#     E(g)$label <- W
# 
#     gr <- igraph_to_networkD3(g)
#     gr$nodes$label <- unlist(labels, use.names=FALSE)
# 
#     labs <- seq(1,length(gr$nodes$label)) %>% map_chr(~as.character(.))
#     if(!is.null(original_dataset)){
#         matching_samples <- get_samples_lables_matches(gr, original_dataset)
#         for(i in seq(1,length(gr$nodes$label))){
#             if(length(strsplit(matching_samples[i], ",")[[1]]) <= 4){
#                 labs[i] <- matching_samples[i]
#             }
#         }
#     }
#     
#     if(is.null(original_dataset)){
#         matching_samples <- ""
#     }else{
#         matching_samples <- paste0("<br></p></font><p>Matching:<br>", matching_samples)
#     }
#    
#     nodes <- data.frame(
#         id = gr$nodes$name,
#         label = labs, #gr$nodes$label,
#         shadow = TRUE,
#         shape = "box",
#         title = paste0("<font color=\"black\"><p>Node:<br>",
#                     gr$nodes$label,
#                     matching_samples)
#     )
#     # compute genotype difference from source to destination
#     diff <- map2(gr$links$source, gr$links$target,
#                 function(x, y){
#                     a <- unlist(strsplit(gr$nodes$label[x+1], ", "))
#                     b <- unlist(strsplit(gr$nodes$label[y+1], ", "))
#                     d <- setdiff(b,a)
#                     d
#                 })
# 
#     edges <- data.frame(
#         from = gr$links$source+1,
#         to = gr$links$target+1,
#         arrows = "to",
#         shadow = TRUE,
#         label = signif(gr$links$value, digits = 4),
#         title = paste0("<font color=\"black\"><p>Edge:<br> from: ",
#                     gr$nodes$label[gr$links$source+1] ,"<br>to : ",
#                     gr$nodes$label[gr$links$target+1] ,
#                     "<br>New mutations:<br>",
#                     diff, "</p></font>")
#     )
# 
#     grp <- toVisNetworkData(g, idToLabel = TRUE)
# 
#     visNetwork(nodes, edges)
# }

#' VisNetwork graph output (default)
#'
#' Draws the output graph using VisNetwork
#'
#' @param out the output object of CIMICE (es, from quick run)
#' @param ... other arguments for format_labels
#'
#' @return visNetwork object representing g as described
#'
#' @examples
#' draw_visNetwork(quick_run(example_dataset()))
#'
#' @export draw_visNetwork
draw_visNetwork <- function(out, ...){
    g <- out$topology 
    W <- out$weights
    labels <- make_labels(out, ...)
    
    V(g)$label <- out$labels
    E(g)$label <- W
    
    gr <- igraph_to_networkD3(g)
    gr$nodes$label <- unlist(out$labels, use.names=FALSE)
    
    nodes <- data.frame(
        id = gr$nodes$name,
        label = labels, #gr$nodes$label,
        shadow = TRUE,
        shape = "box",
        title = paste0("<font color=\"black\"><p>Genotype:<br>",
                       gr$nodes$label, 
                       "<br><p>Matching samples:<br>",
                       out$matching_samples)
    )
    # compute genotype difference from source to destination
    diff <- map2(gr$links$source, gr$links$target,
                 function(x, y){
                     a <- unlist(strsplit(gr$nodes$label[x+1], ", "))
                     b <- unlist(strsplit(gr$nodes$label[y+1], ", "))
                     d <- setdiff(b,a)
                     d
                 })
    
    edges <- data.frame(
        from = gr$links$source+1,
        to = gr$links$target+1,
        arrows = "to",
        shadow = TRUE,
        label = signif(gr$links$value, digits = 4),
        title = paste0("<font color=\"black\"><p>Edge:<br> from:",
                       gr$nodes$label[gr$links$source+1] ,"<br>to : ",
                       gr$nodes$label[gr$links$target+1] ,
                       "<br>New mutations:<br>",
                       diff, "</p></font>")
    )
    
    grp <- toVisNetworkData(g, idToLabel = TRUE)
    
    visNetwork(nodes, edges) # %>% visHierarchicalLayout(sortMethod="directed")
}


# === DEPRECATED ===
#  
# Get associations between graph nodes and the samples of the original dataset
# 
# @param gr a networkD3 graph
# @param original_dataset a dataset with row and column names (samples, genes)
# 
# @return a character vector of comma separated sample identifiers. 
# These lists consists of the identifiers of the samples that match a certain node in the graph.
# (follows node's orderd in the graph).
#  
# @examples
# require(dplyr)
# require(networkD3)
# out <- quick_run(example_dataset())
# gr <- igraph_to_networkD3(out$topology)
# gr$nodes$label <- unlist(out$labels, use.names=FALSE) 
# get_samples_lables_matches(gr, example_dataset())
#
# @export get_samples_lables_matches
# get_samples_lables_matches <- function(gr, original_dataset){
#     # prepare a matrix that will report for each pair (graph_node, sample) if the sample belongs to that graph_node 
#     match_matrix <- matrix(rep(0, times=length(gr$nodes$name)*nrow(original_dataset)), ncol = nrow(original_dataset))
#     # extract the list of considered genes in this analysis
#     considered_genes <- strsplit( paste(gr$nodes$label, collapse = ", "), ",\\s*")[[1]] %>% unique 
#     considered_genes <- considered_genes[considered_genes!="Clonal"]
#     # find pairs (node, sample) with complete matching on wt/mutated genes
#     for(node_i in seq(1,length(gr$nodes$name))){
#         for(sample in seq(1,nrow(original_dataset))){
#             matching <- TRUE
#             present_genes <- strsplit(gr$nodes$label[[node_i]] , ",\\s*")[[1]]
#             for(gene in considered_genes){
#                 if(original_dataset[sample,gene] != gene %in% present_genes){
#                     matching <- FALSE
#                 }
#             }
#             match_matrix[node_i, sample] <- matching
#         }
#     }
#     # prepare strings
#     matching_samples <- rep("", times = length(gr$nodes$name))
#     for(node_i in seq(1,length(gr$nodes$name))){
#         matching_samples[node_i] <- rownames(original_dataset)[which(match_matrix[node_i, ] == TRUE)] %>% paste(collapse = ", ")
#     }
#     matching_samples
# }


#' Dot graph output
#'
#' Represents this graph in dot format (a textual output format)
#'
#' @param out the output object of CIMICE (es, from quick run)
#' @param ... other arguments for format_labels
#'
#' @return a string representing the graph in dot format
#'
#' @examples
#' to_dot(quick_run(example_dataset()))
#'
#' @export to_dot
to_dot <- function(out, ...){
    g <- out$topology 
    W <- out$weights
    labels <- make_labels(out, ...)
    
    V(g)$label <- labels
    E(g)$label <- W

    str <- paste(
        "digraph G{",
        paste(
            imap_chr(V(g)$label, ~ paste( .y, '[label="', .x ,'"]', sep="") ),
            collapse = "\n"),
        paste(
            imap_chr(E(g), ~ paste(tail_of(g, E(g)[[.y]]), " -> ", head_of(g, E(g)[[.y]]), '[label="', round(W[.y], digits = 3) , '"]', sep="" )),
            collapse = "\n"),
        "}",
        sep="\n"
    )

    str
}

