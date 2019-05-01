# Run the VarSelRF function

# Source in the two main function
source('sim_friedman.R')
source('VarSelRF.R')

# Set up a gain function - first one is just rmse
gain_rmse = function(pred, y) sqrt(mean((y - pred)^2))

# Create some data
data_1 = sim_friedman(100)

# Run it through VarSelRF
rf_1 = VarSelRF(X = data_1$X,
                y = data_1$y,
                gain_fun = gain_rmse)
