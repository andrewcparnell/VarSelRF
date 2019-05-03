library(tidyverse)
# Run the VarSelRF function
# Source in the two main function
source('Bruna/VarSelRF/code/sim_friedman.R')
source('Bruna/VarSelRF/code/gain.R')
source('Bruna/VarSelRF/code/predictions.R')
source('Bruna/VarSelRF/code/auxiliar.R')
source('Bruna/VarSelRF/code/VarSelRF.R')

# Create some data
set.seed(123)
data_1 = sim_friedman(100)

X <- cbind(data_1$X, replicate(35, rnorm(n = 100)))
y <- data_1$y
n_trees = 10
n = length(y)
choose_var = floor(sqrt(ncol(X)))
obs_sel = matrix(NA, ncol = n_trees, nrow = n)
feature_sel = matrix(NA, ncol = n_trees, nrow = choose_var)
node_min_size = 10

for(j in 1:n_trees) {
  obs_sel[,j] = sample(1:n, replace = TRUE)
  feature_sel[,j] = sample(1:ncol(X), size = choose_var, replace = FALSE)
}

# Create a list to store all the trees
tree_output = vector(mode = 'list', length = n_trees)
tree_general <- list()
# Using a constant
weights <- rep(0.8, ncol(X))


for(i in 1:n_trees) {
  print(i) 
  tree_output[[i]] = grow_rf_tree(XX = X[obs_sel[,i], feature_sel[,i]], 
                                  yy = y[obs_sel[,i]],
                                  gain_fun = gain_fun,
                                  node_min_size = node_min_size, 
                                  feat = feature_sel[,i], 
                                  weights)
}

#-----------------------------------------------------------------------
tree_output_corr <- list()

fc_corr <- function(var){
  cor(y, var, method = "spearman")
}

weights <- X %>% as.data.frame() %>% 
  purrr::map_dbl(fc_corr) %>% abs()

for(i in 1:n_trees) {
  print(i) 
  tree_output_corr[[i]] = grow_rf_tree(XX = X[obs_sel[,i], feature_sel[,i]], 
                                       yy = y[obs_sel[,i]],
                                       gain_fun = gain_fun,
                                       node_min_size = node_min_size, 
                                       feat = feature_sel[,i], 
                                       weights)
}

#-----------------------------------------------------------------------

