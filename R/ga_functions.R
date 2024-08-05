
# Function for plotting networks
plot_from_matrix <- function(mat, title, file_name){
  output_path <- here('output', 'figures', file_name)
  g <- graph.adjacency(mat, mode = "undirected", weighted = TRUE)
  
  # Set vertex attributes for labels
  V(g)$label <- rownames(mat)
  
  png(output_path, width = 800, height = 600) # Open PNG device
  # Plot the graph with adjusted parameters
  plot(g, layout = layout_with_kk(g),
       vertex.size = 10,            # Increase node size
       vertex.label.cex = 0.8,      # Adjust label size
       vertex.label.dist = 0,     # Increase label distance from nodes
       vertex.color = "lightblue",  # Node color
       vertex.frame.color = "gray", # Node border color
       edge.width = 2,              # Edge thickness
       edge.color = "black",        # Edge color
       
       main = title)
  dev.off()
}

# Function to convert adjacency matrix to edge list
adj_matrix_to_edge_list <- function(adj_matrix) {
  edges <- which(adj_matrix != 0, arr.ind = TRUE)
  # Ensure each edge is represented as (min, max) to handle undirected edges correctly
  edges <- t(apply(edges, 1, function(x) sort(x)))
  # Remove duplicate edges (since undirected)
  edges <- unique(edges)
  return(edges)
}

# Function to find unique edges
find_unique_edges <- function(edges1, edges2) {
  edges1_str <- apply(edges1, 1, paste, collapse = "-")
  edges2_str <- apply(edges2, 1, paste, collapse = "-")
  unique_edges <- setdiff(edges2_str, edges1_str)
  return(unique_edges)
}

# Function to merge adjacency matrices

merge_adjacency_matrices_igraph <- function(matrix1, matrix2) {
  # Convert the adjacency matrices to igraph graph objects
  graph1 <- graph_from_adjacency_matrix(matrix1, mode = "undirected")
  graph2 <- graph_from_adjacency_matrix(matrix2, mode = "undirected")
  
  # Perform the union of the two graphs
  merged_graph <- graph.union(graph1, graph2)
  
  # Convert the merged graph back to an adjacency matrix
  merged_matrix <- as_adjacency_matrix(merged_graph, sparse = FALSE)
  
  return(merged_matrix)
}

# Function to rebuild a matrix from a
reconstruct_symmetric_matrix <- function(upper_entries, n) {
  
  # Create placeholder matrix
  mat <- matrix(0, nrow = n, ncol = n)
  
  # Populate upper entries row by row
  start_idx <- 1
  for (i in 1:(n-1)) {
    end_idx <- start_idx + (n - i - 1)
    mat[i, (i+1):n] <- upper_entries[start_idx:end_idx]
    start_idx <- end_idx + 1
  }
  
  # Fill lower triangle (copy from upper triangle)
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  # Diagonal elements are zero
  diag(mat) <- 0
  
  return(mat)
}

# Function to flatten the upper triangular part of a matrix, rowwise
flatten_upper_triangular_rowwise <- function(mat, n) {
  
  # Initialize an empty vector to store upper triangular elements row-wise
  upper_tri_elements <- numeric(0)
  
  # Loop through each row and column to extract upper triangular elements
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      upper_tri_elements <- c(upper_tri_elements, mat[i, j])
    }
  }
  
  # Return the vector of upper triangular elements row-wise
  return(upper_tri_elements)
}

# Function to fill the the zero bins of a chromosome with the elements of a new chromosome
merge_chromosomes<- function(original_list, new_list){
  
  # Validate input lengths
  if (length(new_list) != sum(original_list == 0)) {
    stop("new_list must be of length equal to the number of zeros in original_list")
  }
  
  # Find indices where original_list is 0
  zero_indices <- which(original_list == 0)
  
  # Replace zeros with elements from new_list
  original_list[zero_indices] <- new_list
  
  return(original_list)
}


# Function to evaluate the fitness of a chromosome
fitness <- function(chromosome){

  # Merging with original fully connected list
  
  bin_vec <- merge_chromosomes(flattened_original_matrix, chromosome)
  
  # Reconstruct the network
  net <- as.GNARnet(reconstruct_symmetric_matrix(bin_vec, 37))
  
  # Fit GNAR on within sample data
  fit <- GNARfit(vts = unemployment_vts[1:275,], net = net, alphaOrder = 1, betaOrder = c(2), globalalpha = TRUE)
  
  # Find next step ahead prediction error
  neg_error <- -1*(sum(unemployment_vts[276,]-predict(fit))^2)
  
  return(neg_error)
  
}

# Function for custom constrained population initialisation

initializePopulation <- function(object, n, ...) {
  k <- 26 # Number of ones
  popSize <- object@popSize
  nBits <- object@nBits
  
  # Ensure k is not greater than nBits
  if (k > nBits) {
    stop("The number of ones (k) cannot be greater than the number of bits (nBits).")
  }
  
  population <- matrix(0, nrow = popSize, ncol = nBits) # Matrix to hold the binary vectors
  
  # For each member of the initial population, k bits are randomly sampled and set equal to one
  for (i in 1:popSize) {
    ones_positions <- sample(1:nBits, k)
    population[i, ones_positions] <- 1
  }
  
  return(population)
  
}

# Function for custom constrained crossover

customCrossover <- function(object, parents, ...) {
  k <- 26
  nBits <- object@nBits
  parents <- object@population[parents, , drop = FALSE] 
  children <- matrix(0, nrow = 2, ncol = nBits) # Matrix to store children
  fitnessChildren <- rep(NA, 2)
  
  crossover_point <- sample(1:(nBits-1), 1) # Crossover point is randomly selected 
  
  # Two offspring are created from the parents
  children[1,] <- c(parents[1,1:crossover_point], parents[2,(crossover_point+1):nBits])
  children[2,] <- c(parents[2,1:crossover_point], parents[1,(crossover_point+1):nBits])
  
  # If the generated chromosome doesnt have k edges, here we randomly add or remove edges to bring the total to k
  
  for (i in 1:2) {
    num_ones <- sum(children[i,])
    if (num_ones > k) {
      # Remove excess ones
      ones_positions <- which(children[i,] == 1)
      remove_positions <- sample(ones_positions, num_ones - k)
      children[i, remove_positions] <- 0
    } else if (num_ones < k) {
      # Add missing ones
      zero_positions <- which(children[i,] == 0)
      add_positions <- sample(zero_positions, k - num_ones)
      children[i, add_positions] <- 1 
      
    }
  }
  
  list(children = children, fitness = fitnessChildren)
}

# Function for custom constrained mutation

customMutation <- function(object, parent, ...) {
  k <- 26
  nBits <- object@nBits
  child <- object@population[parent, , drop = FALSE]
  
  # Perform mutation
  mutation_points <- sample(1:nBits, 2) # Randomly selecting mutation points 
  child[mutation_points] <- 1 - child[mutation_points]
  
  num_ones <- sum(child)
  
  # Correcting to desired number of genes
  
  if (num_ones > k) {
    # Remove excess ones
    ones_positions <- which(child == 1)
    remove_positions <- sample(ones_positions, num_ones - k)
    child[remove_positions] <- 0
  } else if (num_ones < k) {
    # Add missing ones
    zero_positions <- which(child == 0)
    add_positions <- sample(zero_positions, k - num_ones)
    child[add_positions] <- 1
  }
  
  child
}

# Function to create a random initial suggestion

create_row <- function(nCols, nOnes) {
  row <- c(rep(1, nOnes), rep(0, nCols - nOnes))
  sample(row)  # Shuffle the row
}