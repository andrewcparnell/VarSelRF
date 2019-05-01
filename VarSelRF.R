# Function to run Random Forests with a Gain function provided as an argument to allow flexible variable selection

# Main function
VarSelRF = function(X, # feature matrix
                    y, # response variable (vector),
                    gain_fun, # Gain function
                    n_trees = 50, # Number of trees
                    node_min_size = 10, # Minimum size of terminal node to grow
                    choose_var = floor(sqrt(ncol(X))) # Number of features to select for each split
) {
  
  # Algorithm runs as follows
  # 1. Create n_trees data sets boostrapping the rows and with choose_var columns
  # 2. For each data set:
  #   a. Take each variable and each value of X in turn and calculate the gain function
  #   b. Split a variable at the maximum gain
  #   c. Stop when you reach no more possible splits
  #   d. Repeat until tree is maximally grown

  # Step 1 create the indices for each bootstrap data set
  n = length(y)
  obs_sel = matrix(NA, ncol = n_trees, nrow = n)
  feature_sel = matrix(NA, ncol = n_trees, nrow = choose_var)
  for(j in 1:n_trees) {
    obs_sel[,j] = sample(1:n, replace = TRUE)
    feature_sel[,j] = sample(1:ncol(X), size = choose_var, replace = FALSE)
  }
  
  # Create a list to store all the trees
  tree_output = vector(mode = 'list', length = n_trees)
  
  # Set 2 call the function to grow the trees based on the boostrapped data sets 
  # This loop should be parallelised in future
  for(i in 1:n_trees) {
    tree_output[[i]] = grow_rf_tree(X[obs_sel[,i], feature_sel[,i]], 
                                    y[obs_sel[,i]],
                                    gain_fun = gain_fun,
                                    node_min_size = node_min_size)
  }
  
  return(tree_output)
  
}

# Function to grow trees
grow_rf_tree = function(XX, # Bootstrapped features
                        yy, # Bootstrapped response
                        gain_fun,
                        node_min_size) {
  
  # Start by creating a stump
  tree = create_stump(XX, yy)
  
  # Start a while loop that keeps growing trees until the gain is 0
  still_going = TRUE
  while(still_going) {
    # Find all the terminal nodes in the current tree and grow them if the minimum size is big enough
    which_can_grow = which(tree$tree_matrix[,'terminal'] == 1 & tree$tree_matrix[,'node_size'] > node_min_size)
    
    # Always grow the first node that can be grown
    # Find the maximum gain across each of the columns and return the split variable and value associated with it
    best_split = find_max_gain(XX, yy, gain_fun, tree, which_can_grow[1])
    
    # Now split the tree based on that value
    tree = grow_tree(XX, yy, tree, best_split)
    
    # If there are no values that can be grown then stop
    if(length(which_can_grow) == 0) still_going = FALSE
  }
  
  return(tree)

}

# Function to create a stump
create_stump = function(X, y) {
  
  # Each tree has 6 columns and 2 elements
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
  # Gain
  
  # Set up the tree to have two elements in the list as described above
  tree = vector('list', length = 2)
  # Give the elements names
  names(tree) = c('tree_matrix', 
                  'node_indices')
  
  # Create the two elements: first is a matrix
  tree[[1]] = matrix(NA, ncol = 9, nrow = 1)
  
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
    'node_size',
    'gain' # Start this off at 1 so that we can grow the tree
  )
  
  # Set values for stump 
  tree[[1]][1,] = c(1, 1, NA, NA, NA, NA, mean(y), length(y), 1)

return(tree)
  
} # End of function

# Function to find the maximum gain for the current Xs based on a single terminal node
find_max_gain = function(XX,
                         yy,
                         gain_fun,
                         tree,
                         node_to_grow) {

  # First get the predictions from the current tree
  tree_pred = get_predictions(tree, XX)

  # Create somewhere to store the gains for each x value as well
  gain_store = matrix(NA, ncol = 3, nrow = ncol(XX))
  colnames(gain_store) = c('split_var', 'split_val', 'gain_val')

  # Now loop through each of the columns and each of the unique values in each column to find best split
  for(j in 1:ncol(XX)) {
    curr_X = XX[,j]
    unique_X = sort(unique(curr_X))

    # Loop through each unique value of X
    for(i in 1:length(unique_X)) {
      curr_split_var = j
      curr_split_val = unique_X[i]
      # Try growing a tree with this split var/val combination
      try_tree = grow_tree(XX, yy, tree, node_to_grow, curr_split_var, curr_split_val)
      # Get the predictions from this new tree
      try_tree_pred = get_predictions(try_tree)
      # Evaluate the gain
      curr_gain = gain_fun(tree_pred, curr_tree_pred)

      # Fill in the gain_store if it's an improvement
      if(j > 1) {
        if(curr_gain > max(gain_store[,3])) gain_store[j,] = c(j, unique_X[i], curr_gain)
      } else {
        gain_store[1,] = c(1, unique_X[1], 0)
      }
    }

  }

}


# Gets the predicted values from a current set of trees
get_predictions = function(tree, X) {
  
    # If just a single row prediction is easy!  
    if(nrow(tree$tree_matrix) == 1) {
      predictions = rep(tree$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(tree$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(tree, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] = 
          tree$tree_matrix[unique_node_indices[i], 'mu']
      }
    }

  return(predictions)
}

# The fill tree details function takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in ecah terminal node
fill_tree_details = function(curr_tree, X) {
  
  # This code should essentially start from ignorance - no indices just a tree
  # Fill in the number of observations and the node indices
  browser()
  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix
  
  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix
  
  # Start with dummy node indices
  node_indices = rep(1, nrow(X))
  
  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {
    # Get the parent
    curr_parent = tree_matrix[i,'parent']
    
    # Find the split variable and value of the parent
    split_var = tree_matrix[curr_parent,'split_variable']
    split_val = tree_matrix[curr_parent, 'split_value']
    
    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table
  
  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))  
  
} # End of function

grow_tree = function(XX, yy, curr_tree, node_to_split, split_variable, split_value) {
  
  browser()
  # Set up holder for new tree
  new_tree = curr_tree
  
  # Get the list of terminal nodes
  terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1) # Create the list of terminal nodes
  
  # Find terminal node sizes
  terminal_node_size = new_tree$tree_matrix[terminal_nodes,'node_size']
  
  # Add two extra rows to the tree in question
  new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                               c(1, NA, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                               c(1, NA, NA, NA, NA, NA, NA, NA, NA))
  
  curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
  new_tree$tree_matrix[node_to_split,1:6] = c(0, # Now not temrinal
                                              nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                              nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                              curr_parent,
                                              split_variable,
                                              split_value)
  
  #  Fill in the parents of these two nodes
  new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split 
  new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split 
  
  # Now call the fill function on this tree
  new_tree = fill_tree_details(new_tree, X)
  
  # Return new_tree
  return(new_tree)
  
} # End of grow_tree function
