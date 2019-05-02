# Run the VarSelRF function

# Source in the two main function
source('sim_friedman.R')
source('VarSelRF.R')

# Call in foreach library
library(foreach)
library(doSNOW)
library(randomForest)

# Set up a gain function - first one is just reduction in rmse
rmse = function(pred, y) sqrt(mean((y - pred)^2))
gain_rmse = function(new_pred, old_pred, y) sqrt(mean((y - old_pred)^2)) - sqrt(mean((y - new_pred)^2))

# Create some data
set.seed(123)
data_1 = sim_friedman(1000, p = 20)

# Run it through VarSelRF
rf_1 = VarSelRF(X = data_1$X,
                y = data_1$y,
                gain_fun = gain_rmse,
                n_trees = 20)

# Try it against some new data
data_new = sim_friedman(500, p = 0)

# Get predictions from the model
pred_new = predict_VarSelRF(rf_1, newX = data_new$X)

# Quick plot
plot_range = range(c(data_new$y, pred_new))
plot(data_new$y, pred_new, xlim = plot_range, ylim = plot_range)
abline(a=0, b=1)

# # Set up a gain function - first one is just reduction in rmse
# gain_with_depth = function(new_pred, old_pred, y, tree) {
#   depth_of_current_split = nrow(tree$tree_matrix) 
#   out = (sqrt(mean((y - new_pred)^2)) - sqrt(mean((y - old_pred)^2)))^depth_of_current_split
# }
#   
# 

# Compare with randomForest
rf_2 = randomForest(x = data_1$X,
                    y = data_1$y,
                    ntree = 20)

# Look at a tree
getTree(rf_2, 1)

# Get predictions from the model
pred_new_2 = predict(rf_2, newdata = data_new$X)

# Quick plot
plot_range = range(c(data_new$y, pred_new_2))
plot(data_new$y, pred_new_2, xlim = plot_range, ylim = plot_range)
abline(a=0, b=1)
