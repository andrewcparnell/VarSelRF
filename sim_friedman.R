
# Simulate some data according to the Friedman example from the MARS paper
sim_friedman = function(n, p = 0, res_sd = 0.5, pars = c(10, 20, 10, 5)) {
  # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
  X = matrix(NA, nrow = n, ncol = 5 + p)
  for(i in 1:ncol(X)) X[,i] = rnorm(n, 0, 1)
  y = mean = rep(NA, n)
  err = rnorm(n, sd = res_sd)
  mean = pars[1]*sin(pi*X[,1]*X[,2]) + pars[2] * (X[,3]-0.5)^2 + pars[3] * X[,4] + pars[4] * X[,5]
  y = mean + err
  return(list(y = y, X = X, mean = mean, res_sd = res_sd))
}
