# Function to run Random Forests with a Gain function provided as an argument to allow flexible variable selection

#' reg_rf
#'
#' Performs regularisation in random forests
#'
#' @param X The feature matrix, 
#' @param y The response variables,
#' @param gain_fun The gain function used, 
#' @param n_trees The number of trees in the RF, 
#' @param node_min_size The minimum size of terminal node to grow, 
#' @param choose_var  The number of features to select for each split. 
#' Algorithm runs as follows
#' 1. Create n_trees data sets boostrapping the rows and with choose_var columns
#' 2. For each data set:
#'  a. Take each variable and each value of X in turn and calculate the gain function
#'   b. Split a variable at the maximum gain
#'  c. Stop when you reach no more possible splits
#'  d. Repeat until tree is maximally grown
#'
#' @return The resulting model
#'
#' @export
#' 

# reg_rf <-  function(X, # feature matrix
#                     y, # response variable (vector),
#                     gain_fun, # Gain function
#                     n_trees = 50, # Number of trees
#                     node_min_size = 10, # Minimum size of terminal node to grow
#                     choose_var = floor(sqrt(ncol(X))), # Number of features to select for each split
#                     ... ) {
#   
#   # Step 1 create the indices for each bootstrap data set
#   n = length(y)
#   obs_sel = matrix(NA, ncol = n_trees, nrow = n)
#   feature_sel = matrix(NA, ncol = n_trees, nrow = choose_var)
#   
#   for(j in 1:n_trees) {
#     obs_sel[,j] = sample(1:n, replace = TRUE)
#     feature_sel[,j] = sample(1:ncol(X), size = choose_var, replace = FALSE)
#   }
#   
#   # Create a list to store all the trees
#   tree_output = vector(mode = 'list', length = n_trees)
#   
#   # Set 2 call the function to grow the trees based on the boostrapped data sets 
#   # This loop should be parallelised in future
#   # Also should put a progress bar in here 
#   for(i in 1:n_trees) {
#     print(i) 
#     tree_output[[i]] = grow_rf_tree(XX = X[obs_sel[,i], feature_sel[,i]], 
#                                     yy = y[obs_sel[,i]],
#                                     gain_fun = gain_fun,
#                                     node_min_size = node_min_size, 
#                                     feat = feature_sel[,i], 
#                                     weights)
#   }
#   
#   return(list(trees = tree_output,
#               y = y,
#               X = X,
#               n_trees = n_trees,
#               gain_fun = gain_fun,
#               node_min_size = node_min_size))
#   
# }

#' 
#' grow_rf_tree
#'
#' Function to grow trees
#'
#' @param XX The bootstrapped matrix of predictors,
#' @param yy The bootstrapped response,
#' @param gain_fun The gain function used, 
#' @param node_min_size The minimum size of terminal node to grow, 
#'
#' @return The resulting model
#'
#' @export
#' 
grow_rf_tree <-  function(XX, # Bootstrapped features
                          yy, # Bootstrapped response
                          gain_fun,
                          node_min_size,
                          feat,
                          ...) {
  
  # Start by creating a stump
  tree = create_stump(XX, yy)
  used_vars <- gains <-  vector(length = 0)
  # Start a while loop that keeps growing trees until the gain is 0
  still_going = TRUE
  #count = 1
  while(still_going) {
    # Find all the terminal nodes in the current tree and grow them 
    # if the minimum size is big enough
    which_can_grow = which(
      tree$tree_matrix[,'terminal'] == 1 & 
        tree$tree_matrix[,'node_size'] > node_min_size)
    
    if(length(which_can_grow) == 0) {
      # If there are no values that can be grown then stop
      still_going = FALSE
    } else { 
      # Always grow the first node that can be grown
      # Find the maximum gain across each of the columns and return 
      # the split variable and value associated with it
      best_split = find_max_gain(XX, yy, 
                                 gain_fun, 
                                 tree, 
                                 node_to_grow = which_can_grow[1],
                                 current_used_vars = used_vars,
                                 feat = feat, 
                                 weights)
      
      used_vars <- c(used_vars, best_split[4])
      gains <- c(gains, best_split[3])
      # Now split the tree based on that value
      tree = grow_tree(XX, yy, 
                       curr_tree = tree, 
                       node_to_split = which_can_grow[1],
                       split_variable = best_split[1],
                       split_value = best_split[2])
    }
  }
  
  return(list(tree = tree, 
              gains = gains, 
              used_vars = used_vars))
}


#' find_max_gain
#'
#' Function to find the maximum gain for the current Xs based on a single terminal node
#'
#' @param XX The bootstrapped matrix of predictors,
#' @param yy The bootstrapped response,
#' @param gain_fun The gain function used, 
#' @param tree The current tree,
#' @param node_to_grow The node to grow from, 
#' @param weights The weights for the variables
#'
#' @return 
#'
#' @export
#' 
#' 
find_max_gain <-  function(XX,
                           yy,
                           gain_fun,
                           tree,
                           node_to_grow, 
                           weights, 
                           current_used_vars,
                           feat,
                           ...) {
  XX <- as.matrix(XX)
  
  # First get the predictions from the current tree
  tree_pred = get_predictions(tree, XX)
  
  # Create somewhere to store the gains for each x value as well
  gain_store = matrix(NA, ncol = 4, nrow = ncol(XX))
  colnames(gain_store) = c('split_var', 'split_val', 'gain_val', 'var_sel')
  
  # Now loop through each of the columns and each of the unique values
  # in each column to find best split  
  for(j in 1:ncol(XX)){
    curr_X = XX[,j]
    var_sel = feat[j]
    # Remember the only X values you can split on are those 
    # in that terminal node
    unique_X = sort(unique(curr_X[tree$node_indices == node_to_grow]))
    
    # If there's only one X value stop right now
    if(length(unique_X) == 1) {
      gain_store[j,] = c(1, unique_X, 0, 0)
    } else {   
      # Loop through each unique value of X
      for(i in 2:length(unique_X)) {
        curr_split_var = j 
        curr_split_val = unique_X[i]
        # Try growing a tree with this split var/val combination
        try_tree = grow_tree(XX, yy, tree, node_to_grow, 
                             curr_split_var, curr_split_val)
        # Get the predictions from this new tree
        try_tree_pred = get_predictions(try_tree, XX)
        # Evaluate the gain
        curr_gain = gain_fun(new_pred = tree_pred, 
                             old_pred = try_tree_pred, 
                             y = yy, 
                             var_sel_ind = var_sel, 
                             current_used_vars = current_used_vars, 
                             weights)
        
        # Fill in the gain_store if it's an improvement
        if(i > 2) {
          if(curr_gain > gain_store[j,3])
            { gain_store[j,] = c(j, unique_X[i], curr_gain, var_sel) } 
        } else {
          gain_store[j,] = c(j, curr_split_val, curr_gain, var_sel)
        }
      } # End of loop through unique X values
    } # End of if statment on length of unique X
  } # End of loop through features
  
  # Find the best gain and return that row of it
  best_gain <- gain_store[which.max(gain_store[,3]),]
  return(best_gain)
  
} # End of function


