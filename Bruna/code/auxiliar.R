#' create_stump
#'
#' Function to create a stump
#'
#' @param X The matrix of predictors for the stump,
#' @param y The response variable for the stump,
#'
#' @return A stump to be filled with the model 
#'
#' @export
#' 

create_stump <-  function(X, y) {
  
  # Each tree has 8 columns and 2 elements
  # The 3 elements are the tree matrix, and the node indices
  # The tree matrix has columns:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu values
  # Node size
  
  # Set up the tree to have two elements in the list as described above
  tree = vector('list', length = 2)
  # Give the elements names
  names(tree) = c('tree_matrix', 
                  'node_indices')
  
  # Create the two elements: first is a matrix
  tree[[1]] = matrix(NA, ncol = 8, nrow = 1)
  
  # Second is the assignment to node indices
  tree[[2]] = rep(1, length(y))
  
  # Create column names
  colnames(tree[[1]]) = c(
    'terminal',
    'child_left',
    'child_right',
    'parent',
    'split_variable',
    'split_value',
    'mu',
    'node_size'
  )
  
  # Set values for stump 
  tree[[1]][1,] = c(1, 1, NA, NA, NA, NA, mean(y), length(y))
  
  return(tree)
  
} 

#' fill_tree_details
#'
#' Takes a tree matrix and returns the number of obs in each node in
#'  it and the indices of each observation in each terminal node
#'
#' @param curr_tree The current tree,
#' @param X The matrix of predictors,
#'
#' @return 
#'
#' @export
#' 
#' 
fill_tree_details <- function(curr_tree, X) {
  
  # This code should essentially start from ignorance - no indices,
  # just a tree
  # Fill in the number of observations and the node indices
  
  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix
  
  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix
  X <- as.matrix(X)
  # Start with dummy node indices
  node_indices = rep(1, nrow(X))
  
  # For all but the top row, find the number of observations 
  # falling into each one 
  for(i in 2:nrow(tree_matrix)) {
    # Get the parent
    curr_parent = tree_matrix[i, "parent"]
    
    # Find the split variable and value of the parent
    split_var = tree_matrix[curr_parent,'split_variable']
    split_val = tree_matrix[curr_parent, 'split_value']
    
    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == 
                                               curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == 
                                                    curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == 
                                               curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == 
                                                    curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table
  
  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))  
  
} # End of function

grow_tree = function(XX, yy, curr_tree, node_to_split, 
                     split_variable, split_value) {
  
  # Set up holder for new tree
  new_tree = curr_tree
  
  # Get the list of terminal nodes
  terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1) 
  # Create the list of terminal nodes
  
  # Find terminal node sizes
  terminal_node_size = new_tree$tree_matrix[terminal_nodes,'node_size']
  
  # Add two extra rows to the tree in question
  new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                               c(1, NA, NA, NA, NA, NA, NA, NA), 
                               # Make sure they're both terminal
                               c(1, NA, NA, NA, NA, NA, NA, NA))
  
  curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] 
  # Make sure to keep the current parent in there. Will be NA if at the root node
  new_tree$tree_matrix[node_to_split,1:6] = c(0, 
                                              # Now not temrinal
                                              nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                              nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                              curr_parent,
                                              split_variable,
                                              split_value)
  
  #  Fill in the parents of these two nodes
  new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split[1] 
  new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split 
  
  # Now call the fill function on this tree
  new_tree = fill_tree_details(new_tree, XX)
  
  # Finally fill in the mu values for each terminal node
  new_tree = fill_mu(new_tree, yy)
  
  # Fill in the predicted means
  # new_tree_matrix[i,'mu'] = mean(y[node_indices == i])
  # if(any(is.nan(new_tree_matrix[,'mu']))) browser()
  
  
  # Return new_tree
  return(new_tree)
  
} # End of grow_tree function

#' fill_mu
#'
#' Fill in the terminal node values
#'
#' @param tree The tree,
#' @param y The response variables,
#'
#' @return 
#'
#' @export
#' 
fill_mu = function(tree, y) {
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  for(i in which_terminal) {
    tree$tree_matrix[i,'mu'] = mean(y[tree$node_indices == i])
  }
  return(tree)
}
