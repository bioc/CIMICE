utils::globalVariables(c("from", "label", "matching_sample", "rows", "to", "tot", "weight", "weight.x", "weight.y"))

#' Create a stub for a generator 
#'
#' Create a generator topology directly from a dataset.
#' The topology will follow the subset relation.
#' 
#' @param dataset A compacted CIMICE dataset
#'
#' @return a generator, with weight = 0 for all the edges
#' 
#' @examples 
#' make_generator_stub(example_dataset())
#' 
#' @export make_generator_stub
make_generator_stub <- function(dataset){
  graph_data <- example_dataset() %>% dataset_preprocessing()
  generator <- graph_non_transitive_subset_topology(graph_data$samples, graph_data$labels) %>% 
    as_tbl_graph() %>% 
    mutate(matching_sample = graph_data$matching_samples)
  for(label in generator %>% pull(label) ){
    generator <- generator %>% bind_edges(from = label, to = label, node_key = "label")
  }
  
  generator <- generator %E>% mutate(weight = 0)
  generator <- activate(generator, what="nodes")
  list(graph = generator, data = graph_data)
}

#' Add edge weights to a generator
#' 
#' @param generator a generator 
#' @param f_t_w_list a list of triplets (from, to, list), the triplets must not 
#' be nested in the list. For example list("A","B",0.3, "B", "C", 0.2) is a valid
#' input.
#' @param by "labels" or "samples" to use gene labels or sample labels as references for
#' edge identifiers.
#' 
#' 
#' @return the generator with the modified edges (invalid edges are ignored)
#'   
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'  make_generator_stub() %>% 
#'  set_generator_edges(
#'    list(
#'     "D", "A, D", 1 , 
#'     "A", "A, D", 1 , 
#'     "A, D", "A, C, D", 1 , 
#'     "A, D", "A, B, D", 1 , 
#'     "Clonal", "D", 1 , 
#'     "Clonal", "A", 1 , 
#'     "D", "D", 1 , 
#'     "A", "A", 1 , 
#'     "A, D", "A, D", 1 , 
#'     "A, C, D", "A, C, D", 1 , 
#'     "A, B, D", "A, B, D", 1 , 
#'     "Clonal", "Clonal", 1 
#'   ))
#'    
#' @export set_generator_edges      
set_generator_edges <- function(generator, f_t_w_list, by = "labels"){
  
  # create columns for the data frame
  type <- 1
  from_l <- list()
  to_l <- list()
  weight_l <- c()
  for(el in f_t_w_list){
    if(type == 1){
      from_l <- c(from_l, el)
      type <- 2
    }else if(type == 2){
      to_l <- c(to_l, el)
      type <- 3
    }else{
      assert_that(is.numeric(el), msg="weights must be numeric")
      weight_l <- c(weight_l, el)
      type <- 1
    }
  }
  assert_that(length(from_l) == length(to_l) & length(to_l) == length(weight_l), 
              msg="The sequence must be multiple of three (from, to, weight)")
  
  # get node ids by label/sample_label
  if(by == "samples"){
    id_col <- "matching_samples"
  }else{
    id_col <- "label"
  }
  
  from_l <- map_int(from_l, ~ ifelse(is.numeric(.),
                                     ., 
                                     which((generator$graph %N>% pull(!!id_col)) == .)))
  to_l <- map_int(to_l, ~ ifelse(is.numeric(.),
                                     ., 
                                     which((generator$graph %N>% pull(!!id_col)) == .)))
  
  # create the weight dataframe
  weights_df <- data.frame(from=from_l,to=to_l,weight=weight_l)
  
  # add weights
  generator$graph <- generator$graph %E>% left_join(weights_df, by=c("from","to")) %>% 
    mutate(weight = ifelse(is.na(weight.y), weight.x, weight.y)) %>% 
    select(-weight.x, -weight.y)
  
  generator$graph <- activate(generator$graph, what="nodes")
  generator
}

#' Prepare a command to add edge weights to a generator
#' 
#' Prints a string in the form of the command that 
#' sets weights for all the edges of this generator.
#'
#' @param generator a generator 
#' @param by "labels" or "samples" to use gene labels or sample labels as references for
#' edge identifiers.
#' 
#' @return NULL (the string with the function calls is printed on the stdout) 
#'
#' @examples 
#' require(dplyr)
#' example_dataset() %>% 
#'   make_generator_stub() %>%
#'   prepare_generator_edge_set_command()
#'
#' @export prepare_generator_edge_set_command
prepare_generator_edge_set_command <- function(generator, by="labels"){
  edges <- generator$graph %E>% data.frame()
  generator$graph %N>% print
  if(by == "samples"){
    names_vec <- generator$graph %N>% pull(matching_sample)
  }else{
    names_vec <- generator$graph %N>% pull(label)
  }
  out <- c("set_generator_edges(\n  list(\n")
  for(row in seq(1,nrow(edges))){
    out <- c(out, "    ")
    out <- c(out, paste(
      paste('"', names_vec[edges[row, "from"]], '"', sep=""),
      paste('"', names_vec[edges[row, "to"]], '"', sep=""),
      "_", sep=", "))
    if(row < nrow(edges)){
      out <- c(out, ",")
    }
    out <- c(out, "\n")
  }
  out <- c(out, "))\n")
  paste(out) %>% cat
}

#' Plot a generator
#'
#' Simple ggraph interface to draw a generator 
#' 
#' @param generator a generator
#' 
#' @return a basic plot of this generator
#'
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'  make_generator_stub() %>% 
#'  set_generator_edges(
#'    list(
#'     "D", "A, D", 1 , 
#'     "A", "A, D", 1 , 
#'     "A, D", "A, C, D", 1 , 
#'     "A, D", "A, B, D", 1 , 
#'     "Clonal", "D", 1 , 
#'     "Clonal", "A", 1 , 
#'     "D", "D", 1 , 
#'     "A", "A", 1 , 
#'     "A, D", "A, D", 1 , 
#'     "A, C, D", "A, C, D", 1 , 
#'     "A, B, D", "A, B, D", 1 , 
#'     "Clonal", "Clonal", 1 
#'   )) %>% 
#'   finalize_generator %>% 
#'   plot_generator
#'  
#' @export plot_generator
plot_generator <- function(generator){
  ggraph(generator$graph, layout = "sugiyama") +
    geom_node_point() +
    geom_node_label(aes(label=label)) +
    geom_edge_link(
      aes(label=weight, start_cap = label_rect(node1.label),
          end_cap = label_rect(node2.label)),
      arrow = arrow(length = unit(4, 'mm')),
      label_dodge = unit(2.5, 'mm'),
      angle_calc = 'along'
    ) +
    geom_edge_loop(
      aes(label=weight, start_cap = label_rect(node1.label),
          end_cap = label_rect(node2.label)),
      arrow = arrow(length = unit(4, 'mm')),
      label_dodge = unit(2.5, 'mm'),
      angle_calc = 'along'
    ) 
}

#' Finalize generator normalizing edge weights
#' 
#' Checks if a generator can be normalized so that it 
#' actually is a Markov Chain 
#' 
#' @param generator a generator 
#' 
#' @return A generator with edge weights that respect DTMC definition
#' 
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'  make_generator_stub() %>% 
#'  set_generator_edges(
#'    list(
#'     "D", "A, D", 1 , 
#'     "A", "A, D", 1 , 
#'     "A, D", "A, C, D", 1 , 
#'     "A, D", "A, B, D", 1 , 
#'     "Clonal", "D", 1 , 
#'     "Clonal", "A", 1 , 
#'     "D", "D", 1 , 
#'     "A", "A", 1 , 
#'     "A, D", "A, D", 1 , 
#'     "A, C, D", "A, C, D", 1 , 
#'     "A, B, D", "A, B, D", 1 , 
#'     "Clonal", "Clonal", 1 
#'   )) %>% 
#'   finalize_generator
#' 
#' @export finalize_generator
finalize_generator <- function(generator){
  normalization_df <- generator$graph %E>% 
    data.frame() %>% 
    group_by(from) %>% 
    summarise(tot = sum(weight)) 
  
  if(length(which(normalization_df %>% pull(tot) == 0)) != 0){
    stop("Some node has no exiting edge with probability greater than 0.\nCheck your call to set_generator_edges")
  }
  
  generator$graph <- generator$graph %E>% 
    left_join(normalization_df, by="from") %>% 
    mutate(weight = weight/tot) %>% 
    select(-tot)
  
  generator
}

#' Create datasets from generators
#' 
#' Simulate the DTMC associated to the generator to create a dataset
#' that reflects the genotypes of `times` cells, sampled after
#' `time_ticks` passages.
#' 
#' @param generator a generator 
#' @param time_ticks number of steps (updates) of the DTMC associated to the generato
#' @param times number of sumlated cells
#' @param starting_label node from which to start the simulation
#' @param by "labels" or "samples" to use gene labels or sample labels as references to identify the `starting_label`'s node
#' @param mode "full" to generate a matrix with `times` genotypes, "compacted" to *efficiently* create an already compacted dataset 
#' (a dataset showing the genotypes and their respective frequencies)
#'
#' @return the simulated dataset
#'
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'   make_generator_stub() %>% 
#'   set_generator_edges(
#'     list(
#'       "D", "A, D", 1 , 
#'       "A", "A, D", 1 , 
#'       "A, D", "A, C, D", 1 , 
#'       "A, D", "A, B, D", 1 , 
#'       "Clonal", "D", 1 , 
#'       "Clonal", "A", 1 , 
#'       "D", "D", 1 , 
#'       "A", "A", 1 , 
#'       "A, D", "A, D", 1 , 
#'       "A, C, D", "A, C, D", 1 , 
#'       "A, B, D", "A, B, D", 1 , 
#'       "Clonal", "Clonal", 1 
#'   )) %>% 
#'   finalize_generator %>% 
#'   simulate_generator(3, 10)
#'
#' @export simulate_generator
simulate_generator <- function(generator, time_ticks, times, starting_label = "Clonal", by = "labels", mode="full"){
  m <- sparseMatrix(i = generator$graph %E>% pull(from),
               j = generator$graph %E>% pull(to),
               x = generator$graph %E>% pull(weight)) %>% 
    as.matrix()
  # quickly compute n_step matrix (O(log2(times))) matrix multiplications
  m <- m %^% time_ticks
  if(by == "samples"){
    names_vec <- generator$graph %N>% pull(matching_sample)
  }else{
    names_vec <- generator$graph %N>% pull(label)
  }
  v <- m[which(names_vec == starting_label),]
  
  if(mode == "full"){
    # this option easily allows to introduce FP/FN rates and missing data
    dataset <- sample(seq(1,length(v)), times, prob = v, replace = TRUE)
    generator$data$samples[dataset,]
  }else if(mode == "compacted"){
    # this option efficiently computes very large datasets
    dataset <- sample(seq(1,length(v)), times, prob = v, replace = TRUE)
    selected_genotypes <- dataset %>% unique()
    count_df <- data.frame(rows = dataset) %>% count(rows)
    sorted_count_df <- inner_join(data.frame(rows = selected_genotypes), count_df, by="rows") 
    out_matrix <- generator$data$samples[selected_genotypes,]
    list(matrix = out_matrix, 
         counts = sorted_count_df %>% pull(n),
         row_names = rownames(out_matrix))
  }else{
    stop("Invalid mode ", mode, " for simulate_generator")
  }
  
}

#' Perturbate a boolean matrix
#'
#' Given a boolean matrix, randomly 
#' add False Positives (FP), False Negatives (FN) and Missing data 
#' following user defined rates. In the final matrix, missing 
#' data is represented by the value 3.
#' 
#' Note that CIMICE does not support dataset with missing data natively, 
#' so using MIS_rate != 0 will then require some pre-processing.
#'
#' @param dataset a matrix/sparse matrix
#' @param FP_rate False Positive rate
#' @param FN_rate False Negative rate
#' @param MIS_rate Missing Data rate
#' 
#' @return the new, perturbed, matrix
#' 
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'   make_generator_stub() %>% 
#'   set_generator_edges(
#'     list(
#'       "D", "A, D", 1 , 
#'       "A", "A, D", 1 , 
#'       "A, D", "A, C, D", 1 , 
#'       "A, D", "A, B, D", 1 , 
#'       "Clonal", "D", 1 , 
#'       "Clonal", "A", 1 , 
#'       "D", "D", 1 , 
#'       "A", "A", 1 , 
#'       "A, D", "A, D", 1 , 
#'       "A, C, D", "A, C, D", 1 , 
#'       "A, B, D", "A, B, D", 1 , 
#'       "Clonal", "Clonal", 1 
#'   )) %>% 
#'   finalize_generator %>% 
#'   simulate_generator(3, 10) %>% 
#'   perturb_dataset(FP_rate = 0.01, FN_rate = 0.1, MIS_rate = 0.12)
#'
#' @export perturb_dataset
perturb_dataset <- function(dataset, FP_rate = 0, FN_rate = 0, MIS_rate = 0){
  apply(dataset, c(1,2), function(x){
    if(rbernoulli(1, MIS_rate)){
      3
    }else if(x == 0){
      if(rbernoulli(1, FP_rate)){
        1
      }else{
        0
      }
    }else if(x == 1){
      if(rbernoulli(1, FN_rate)){
        0
      }else{
        1
      }
    }
  })
}

