#' get_predictions
#' Gets the predicted values from a current set of trees
#'
#' @param tree The current tree,
#' @param X The matrix of predictors,
#' 
#' @return 
#'
#' @export
#'
get_predictions <-  function(tree, X) {
  X <- as.matrix(X)
  # If just a single row prediction is easy!  
  if(nrow(tree$tree_matrix) == 1) {
    predictions = rep(tree$tree_matrix[1, 'mu'], times = nrow(X))
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

#' get_predictions
#' Gets the predicted values for a new dataset
#'
#' @param model The current tree model,
#' @param newX The new matrix of predictors,
#' 
#' @return 
#'
#' @export
#' 
predict_newdata = function(model, newX) {
  # Create predictions based on a new feature matrix
  
  # Create holder for predicted values
  n_new = nrow(newX)
  y_hat_tree = matrix(NA, ncol = n_new, nrow = length(model))
  
  # Now loop through trees and get predictions
  for (i in 1:length(model)) { 
    # Get current set of trees
    curr_tree = model[[i]]$tree
    
    # Use get_predictions function to get predictions
    y_hat_tree[i,] = get_predictions(curr_tree, 
                                     newX)
  }
  y_hat_mat = apply(y_hat_tree, 2, 'mean')
  
  return(y_hat_mat)
  
} # end of predict function


