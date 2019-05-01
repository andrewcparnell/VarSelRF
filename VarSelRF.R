# Function to run Random Forests with a Gain function provided as an argument to allow flexible variable selection

# Main function
VarSelRF = function(X, # feature matrix
                    y, # response variable (vector),
                    gain_fun, # Gain function
                    n_trees = 50, # Number of trees
                    node_min_size = max(10, 0.05*length(y)), # Minimum size of terminal node to grow
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
  # Also should put a progress bar in here
  
  cl = makeCluster(parallel::detectCores())
  registerDoSNOW(cl)
  pb = txtProgressBar(max = n_trees, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  tree_output = foreach(i = 1:n_trees, 
                  .options.snow = opts,
                   .export = c("grow_rf_tree",
                               "find_max_gain",
                               "get_predictions",
                               "grow_tree",
                               "fill_tree_details",
                               "fill_mu",
                               "create_stump")) %dopar%
    {
      result = grow_rf_tree(X[obs_sel[,i], 
                                   feature_sel[,i]],
                                 y[obs_sel[,i]],
                                 gain_fun = gain_fun,
                                 node_min_size = node_min_size)
      return(result)
      }
  close(pb)
  stopCluster(cl) 

  # for(i in 1:n_trees) {
  #   print(i)
  #   tree_output[[i]] = grow_rf_tree(X[obs_sel[,i], feature_sel[,i]],
  #                                   y[obs_sel[,i]],
  #                                   gain_fun = gain_fun,
  #                                   node_min_size = node_min_size)
  # }
  
  # Now go through for each tree and correct the chosen variables
  for(i in 1:n_trees) {
    curr_features = feature_sel[,i]
    tree_output[[i]]$tree_matrix[,'split_variable'] = 
      curr_features[tree_output[[i]]$tree_matrix[,'split_variable']]
  }
  
  return(list(trees = tree_output,
              y = y,
              X = X,
              n_trees = n_trees,
              gain_fun = gain_fun,
              node_min_size = node_min_size))
  
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
  #count = 1
  while(still_going) {
    # Find all the terminal nodes in the current tree and grow them if the minimum size is big enough
    which_can_grow = which(tree$tree_matrix[,'terminal'] == 1 & tree$tree_matrix[,'node_size'] > node_min_size)
    
    if(length(which_can_grow) == 0) {
      # If there are no values that can be grown then stop
      still_going = FALSE
    } else {
      # Always grow the first node that can be grown
      # Find the maximum gain across each of the columns and return the split variable and value associated with it
      best_split = find_max_gain(XX, yy, gain_fun, tree, which_can_grow[1])
      
      # Now split the tree based on that value
      tree = grow_tree(XX, yy, tree, 
                       node_to_split = which_can_grow[1],
                       split_variable = best_split[1],
                       split_value = best_split[2])
    }
    
  }
  
  return(tree)

}

# Function to create a stump
create_stump = function(X, y) {
  
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

  # Now loop through the columns finding the best split
  for(j in 1:ncol(XX)) {
    curr_X = XX[tree$node_indices == node_to_grow,j]
    # Remember the only X values you can split on are those in that terminal node
    #unique_X = sort(unique(curr_X))
    unique_X = sort(unique(quantile(curr_X,
                        probs = seq(0, 1, length = 22)[-c(1,12)])))
    #unique_X = sort(sample(sort(curr_X)[-1], size = 10))

    # If there's only one X value stop right now
    if(length(unique_X) == 1) {
      gain_store[j,] = c(1, unique_X, 0)
    } else {
      # Loop through each unique value of X
      for(i in 2:length(unique_X)) {
        curr_split_var = j
        curr_split_val = unique_X[i]
        # Try growing a tree with this split var/val combination
        try_tree = grow_tree(XX, yy, tree, node_to_grow, curr_split_var, curr_split_val)
        # Get the predictions from this new tree
        try_tree_pred = get_predictions(try_tree, XX)
        # Evaluate the gain
        curr_gain = gain_fun(try_tree_pred, tree_pred, yy)
        
        # Fill in the gain_store if it's an improvement
        if(i > 2) {
          if(curr_gain > gain_store[j,3]) gain_store[j,] = c(j, unique_X[i], curr_gain)
        } else {
          gain_store[j,] = c(j, curr_split_val, curr_gain)
        }
      } # End of loop through unique X values
    } # End of if statment on length of unique X
  } # End of loop through features

  # Find the best gain and return that row of it
  #if(max(gain_store) < 0) browser()
  return(gain_store[which.max(gain_store[,3]),])
  
} # End of function


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
  
  # Set up holder for new tree
  new_tree = curr_tree
  
  # Get the list of terminal nodes
  terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1) # Create the list of terminal nodes
  
  # Find terminal node sizes
  terminal_node_size = new_tree$tree_matrix[terminal_nodes,'node_size']
  
  # Add two extra rows to the tree in question
  new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                               c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                               c(1, NA, NA, NA, NA, NA, NA, NA))
  
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
  new_tree = fill_tree_details(new_tree, XX)
  
  # Finally fill in the mu values for each terminal node
  new_tree = fill_mu(new_tree, yy)

  # Return new_tree
  return(new_tree)
  
} # End of grow_tree function

# Fill in the terminal node values
fill_mu = function(tree, y) {
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  for(i in which_terminal) {
    tree$tree_matrix[i,'mu'] = mean(y[tree$node_indices == i])
  }
  return(tree)
}

predict_VarSelRF = function(VarSelRF, newX) {
  # Create predictions based on a new feature matrix
  # Note that there is minimal error checking in this - newX needs to be right!
  
  # Create holder for predicted values
  n_new = nrow(newX)
  y_hat_tree = matrix(NA, ncol = n_new, nrow = length(VarSelRF$trees))

  # Now loop through trees and get predictions
  for (i in 1:length(VarSelRF$trees)) {
    # Get current set of trees
    curr_tree = VarSelRF$trees[[i]]
    # Use get_predictions function to get predictions
    y_hat_tree[i,] = get_predictions(curr_tree, 
                                     newX)
  }
  y_hat_mat = apply(y_hat_tree, 2, 'mean')

  return(y_hat_mat)
  
} # end of predict function

