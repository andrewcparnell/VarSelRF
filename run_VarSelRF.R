# Run the VarSelRF function

# Source in the two main function
source('sim_friedman.R')
source('VarSelRF.R')

# Set up a gain function - first one is just reduction in rmse
gain_rmse = function(new_pred, old_pred, y) sqrt(mean((y - new_pred)^2)) - sqrt(mean((y - old_pred)^2))

# Create some data
set.seed(123)
data_1 = sim_friedman(100)

# Run it through VarSelRF
rf_1 = VarSelRF(X = data_1$X,
                y = data_1$y,
                gain_fun = gain_rmse)

# Get predictions from the model
pred_1 = predict_VarSelRF(rf_1, newX = data_1$X)

# Quick plot
plot(data_1$y, pred_1)

# # Set up a gain function - first one is just reduction in rmse
# gain_with_depth = function(new_pred, old_pred, y, tree) {
#   depth_of_current_split = nrow(tree$tree_matrix) 
#   out = (sqrt(mean((y - new_pred)^2)) - sqrt(mean((y - old_pred)^2)))^depth_of_current_split
# }
#   
# 
