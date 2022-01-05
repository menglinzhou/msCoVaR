####################################################################
#### Simulation study ##############################################
path = ".../"
source(paste(path, "Functions_CoVaR.R", sep = ""))

#### Take logistic distribution for example ########################
#### For other models, only change truepar, n, and dist ################
level <- c(0.05,0.05)

#### Data generation ###############################################
dist <- "log" 
truepar <- 0.6 ## true parameter
n = 2000  ## sample size
N = 100   ## number of replications
data <- data_generate(real = truepar, n, N, dist = "log")

#### M estimation ##################################################
k_m <- 180 ## selected fllowing Einmahl et al. [2012], aiming at minimizing RMSE
estpar <- matrix(0, nrow = N, ncol = length(truepar))
for(i in 1:N){
  D <- matrix(data[,i], ncol = 2)
  estpar[i] <- M_estimate(X = D, family = dist, m = k_m)
}

#### Sample fraction calculation and tail index estimation #########
k <- gamma_hat <- rep(0, N)
for(i in 1:N){
  D <- matrix(data[,i], ncol = 2)
  re <- index_estimate(dat = D[,2]) ## only for Y
  k[i] <- re[1]
  gamma_hat[i] <- re[2]
}

#### VaR_Y(p) & eta_star estimation #################################
VaR_Y_hat <- eta_star_hat <- rep(0, N)
for(i in 1:N){
  D <- matrix(data[,i], ncol = 2)
  VaR_Y_hat[i] <- VaR_estimate(dat = D[,2], k = rep(k[i],2), p = level[1])
  eta_star_hat[i] <- eta_estimate(par_hat = estpar[i,], p = level, family = dist)
}

#### CoVaR_{Y|X} estimation and real values calculation ##############
CoVaR_hat <- VaR_Y_hat*(eta_star_hat)^(-gamma_hat)

real_values <- real_compute(par = truepar, family = dist, p = level)
