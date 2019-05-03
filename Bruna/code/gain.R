gain_rmse = function(new_pred, old_pred, y, var_sel_ind, 
                     current_used_vars,...){
  gain <- sqrt(mean((y - new_pred)^2)) - sqrt(mean((y - old_pred)^2))
  
  if(!var_sel_ind %in% current_used_vars){ 
    w <- weights[var_sel_ind]
    gain <- w * gain
  }
  return(gain)
} 
gain_fun = gain_rmse
# Create object to save the information about the selected variables 