####################################################################
#### Simulation study ##############################################
path = ".../"
source(paste(path, "functions.R", sep = ""))

level = c(0.05,0.05)  ## (p1, p2)
n = 3000  ## sample size
N = 1000   ## number of replications

#### Logistic distribution ########################

#### Data generation ###############################################
dist = "log" 
truepar = 0.6 ## true parameter
data = data_generate(real = truepar, n, N, dist = dist)

#### M estimation ##################################################
m = 270 ## selected fllowing Einmahl et al. [2012], aiming at minimizing RMSE
estpar = matrix(0, nrow = N, ncol = length(truepar))
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  estpar[i,] = M_estimate(X = D, family = dist, m = m)
}

#### Tail index estimation #########
k1 = 360 
gamma_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  gamma_hat[i] = tail_estimate(dat = D[,2], k=k1) ## only for Y
}

#### VaR_Y(p) & eta_star estimation #################################
k2 = 360
VaR_Y_hat = eta_star_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  VaR_Y_hat[i] = Qvar_esti(dat = D[,2], gamma = gamma_hat[i], 
                           k = k2, p = level[2])
  eta_star_hat[i] = eta_estimate(par_hat = estpar[i,], p = level, family = dist)
}

#### CoVaR_{Y|X} estimation and real values calculation ##############
CoVaR_hat = VaR_Y_hat*(eta_star_hat)^(-gamma_hat)

real_values = real_compute(par = truepar, family = dist, p = level)



#### HR distribution ########################

#### Data generation ###############################################
dist = "hr" 
truepar = 2.5 ## true parameter
data = data_generate(real = truepar, n, N, dist = dist)

#### M estimation ##################################################
m = 420
estpar = matrix(0, nrow = N, ncol = length(truepar))
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  estpar[i,] = M_estimate(X = D, family = dist, m = m)
}

#### Tail index estimation #########
k1 = 410
gamma_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  gamma_hat[i] = tail_estimate(dat = D[,2], k=k1) ## only for Y
}

#### VaR_Y(p) & eta_star estimation #################################
k2 = 420
VaR_Y_hat = eta_star_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  VaR_Y_hat[i] = Qvar_esti(dat = D[,2], gamma = gamma_hat[i], 
                           k = k2, p = level[2])
  eta_star_hat[i] = eta_estimate(par_hat = estpar[i,], p = level, family = dist)
}

#### CoVaR_{Y|X} estimation and real values calculation ##############
CoVaR_hat = VaR_Y_hat*(eta_star_hat)^(-gamma_hat)

real_values = real_compute(par = truepar, family = dist, p = level)


#### Asymmetric logistic distribution ########################

#### Data generation ###############################################
dist = "alog" 
truepar = c(0.6, 0.5, 0.8) ## true parameter
data = data_generate(real = truepar, n, N, dist = dist)

#### M estimation ##################################################
m = 240
estpar = matrix(0, nrow = N, ncol = length(truepar))
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  estpar[i,] = M_estimate(X = D, family = dist, m = m)
}

#### Tail index estimation #########
k1 = 410
gamma_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  gamma_hat[i] = tail_estimate(dat = D[,2], k=k1) ## only for Y
}

#### VaR_Y(p) & eta_star estimation #################################
k2 = 410
VaR_Y_hat = eta_star_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  VaR_Y_hat[i] = Qvar_esti(dat = D[,2], gamma = gamma_hat[i], 
                           k = k2, p = level[2])
  eta_star_hat[i] = eta_estimate(par_hat = estpar[i,], p = level, family = dist)
}

#### CoVaR_{Y|X} estimation and real values calculation ##############
CoVaR_hat = VaR_Y_hat*(eta_star_hat)^(-gamma_hat)

real_values = real_compute(par = truepar, family = dist, p = level)


#### t distribution ########################

#### Data generation ###############################################
dist = "t" 
truepar = c(3, 0.6) ## true parameter
data = data_generate(real = truepar, n, N, dist = dist)

#### M estimation ##################################################
m = 90
estpar = matrix(0, nrow = N, ncol = length(truepar))
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  estpar[i,] = M_estimate(X = D, family = dist, m = m)
}

#### Tail index estimation #########
k1 = 30
gamma_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  gamma_hat[i] = tail_estimate(dat = D[,2], k=k1) ## only for Y
}

#### VaR_Y(p) & eta_star estimation #################################
k2 = 150
VaR_Y_hat = eta_star_hat = rep(0, N)
for(i in 1:N){
  D = matrix(data[,i], ncol = 2)
  VaR_Y_hat[i] = Qvar_esti(dat = D[,2], gamma = gamma_hat[i], 
                           k = k2, p = level[2])
  eta_star_hat[i] = eta_estimate(par_hat = estpar[i,], p = level, family = dist)
}

#### CoVaR_{Y|X} estimation and real values calculation ##############
CoVaR_hat = VaR_Y_hat*(eta_star_hat)^(-gamma_hat)

real_values = real_compute(par = truepar, family = dist, p = level)

