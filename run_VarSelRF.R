# Run the VarSelRF function

# Source in the two main function
source('sim_friedman.R')
source('VarSelRF.R')

# Set up a gain function - first one is just reduction in rmse
gain_rmse = function(pred1, pred2, y) sqrt(mean((y - pred1)^2)) - sqrt(mean((y - pred2)^2))

# Create some data
set.seed(123)
data_1 = sim_friedman(100)

# Run it through VarSelRF
rf_1 = VarSelRF(X = data_1$X,
                y = data_1$y,
                gain_fun = gain_rmse)

# Get predictions from the model
pred_1 = predict_VarSelRF(rf_1, newX = data_1$X)
