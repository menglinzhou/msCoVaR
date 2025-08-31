##############################################################
#### Load Packages ##################################
library(mvtnorm)
library(evd)
library(cubature)
library(parallel)
library(doParallel)
library(foreach)
library(rugarch)
library(evir)
library(sn)
library(MASS)
library(Matrix)
##############################################################


##############################################################
#### Parametric upper tail dependence function
## This function gives the upper TD function for specific dist
## Inputs: X: a vector of bivariate data
##         par: parameters for distribution
##         family: name of distribution: "log" is logistic
##                                       "HR" is Husler-Reiss distribution
##                                       "alog" is asymmetric logistic
##                                       "t" stands for t distirbution
################################################################
tail_dependence = function(X, par, family){
  
  x = X[1]
  y = X[2]
  
  if(family == "log"){
    theta = par[1]
    return(x+y-(x^(1/theta)
                +y^(1/theta))^(theta))}
  
  if(family == "alog"){
    theta = par[1]
    psi_1 = par[2]
    psi_2 = par[3]
    return(psi_1*x+psi_2*y-((x*psi_1)^(1/theta)
                            +(psi_2*y)^(1/theta))^(theta))}
  
  if(family == "t"){
    nu = par[1]
    rho = par[2]
    y1 = sqrt((nu+1)/(1-rho^2))*(rho-(y/x)^(-1/nu))
    y2 = sqrt((nu+1)/(1-rho^2))*(rho-(x/y)^(-1/nu))
    return(x*pt(y1,nu+1)+y*pt(y2,nu+1))}
  
  if(family == "hr"){
    theta = par[1]
    z1 = 1/theta + theta/2*(log(x/y))
    z2 = 1/theta + theta/2*(log(y/x))
    l = x*pnorm(z1) + y*pnorm(z2)
    return(x+y-l)}
}


##############################################################
#### g(x,y) function in M-estimation
## This function generates g(x,y) for specific dist
## Inputs: family: name of distribution: "log" is logistic
##                                       "HR" is Husler-Reiss distribution
##                                       "alog" is asymmetric logistic
##                                       "t" stands for t distirbution
################################################################
generate_gfun = function(family){
  if(family == "alog"){
    
    gfun = function(x){
      return(c(1,x[1], 2*(x[1]+x[2])))}}
  
  if(family == "t"){
    
    gfun = function(x){
      return(c(x[1], (x[1]+x[2])))}
  }
  
  if(family == "log"){
    
    gfun = function(x){1}
    
  }
  
  if(family == "hr"){
    
    gfun = function(x){x[1]}
    
  }
  return(gfun)
}

################################################################
#### Nonparametric upper tail dependence function
## This function gives the nonparametric upper TD function
## Inputs: X: a vector of bivariate data
##         Rank: rank matrix of the dataset
##         m: sample fraction in the nonparametric estimation formula
################################################################
nonpar_tail = function(X, Rank, m){
  ## X stands for the bivariate data
  ## m is the parameter in nonparametric estimation formula
  
  x = nrow(Rank)+0.5-m*X[1]
  y = nrow(Rank)+0.5-m*X[2]
  
  num = sum(ifelse(Rank[,1]>x,1,0)*ifelse(Rank[,2]>y,1,0))

  return(num/m)
}

################################################################
#### Optimization function for M-estimation
## This function gives the nonparametric upper TD function
## Inputs: par: parameters for distribution
##         family: name of distribution: "log" is logistic
##                                       "HR" is Husler-Reiss distribution
##                                       "alog" is asymmetric logistic
##                                       "t" stands for t distirbution
##         gfun: generated function g(x,y)
##         m: sample fraction in the nonparametric estimation formula
##         Rank: rank matrix of the dataset
################################################################
optim_fun = function(par, family, gfun, m, Rank){
  
  gfun = match.fun(gfun)
  
  if(family == "log" || family == "hr"){
    ## one dimension
    fun_0 = function(x) {(tail_dependence(X = x,par[1], family) - nonpar_tail(X = x, Rank, m))*gfun(x)}
    d_0 = cubature::adaptIntegrate(fun_0, rep(0, 2), rep(1, 2), maxEval=10000, tol = 1e-10)$integral
    return(d_0^2)
  }
  
  if(family == "t"){
    ## two dimensions
    fun_0 = function(x) (tail_dependence(X = x,par, family) - nonpar_tail(X = x, Rank, m))*(gfun(x)[1])
    fun_1 = function(x) (tail_dependence(X = x,par, family) - nonpar_tail(X = x, Rank, m))*(gfun(x)[2])
    d_0 = cubature::adaptIntegrate(fun_0, rep(0, 2), rep(1, 2), maxEval=10000, tol = 1e-10)$integral
    d_1 = cubature::adaptIntegrate(fun_1, rep(0, 2), rep(1, 2), maxEval=10000, tol = 1e-10)$integral
    return(d_0^2+d_1^2)
  }
  
  if(family == "alog"){
    ## three dimensions
    fun_0 = function(x) (tail_dependence(X = x,par, family) - nonpar_tail(X = x, Rank, m))*(gfun(x)[1])
    d_0 = cubature::adaptIntegrate(fun_0, rep(0, 2), rep(1, 2), maxEval=10000, tol = 1e-10)$integral
    fun_1 = function(x) (tail_dependence(X = x,par, family) - nonpar_tail(X = x, Rank, m))*(gfun(x)[2])
    d_1 = cubature::adaptIntegrate(fun_1, rep(0, 2), rep(1, 2), maxEval=10000, tol = 1e-10)$integral
    fun_2 = function(x) (tail_dependence(X = x,par, family) - nonpar_tail(X = x, Rank, m))*(gfun(x)[3])
    d_2 = cubature::adaptIntegrate(fun_2, rep(0, 2), rep(1, 2), maxEval=10000, tol = 1e-10)$integral
    return(d_0^2+d_1^2+d_2^2)
  }
}

################################################################
#### Main code for M-estimation
## This function gives M-estimator with a given dist
## Inputs: X: a bivariate dataframe
##         m: sample fraction in the nonparametric estimation formula
##         family: name of distribution: "log" is logistic
##                                       "HR" is Husler-Reiss distribution
##                                       "alog" is asymmetric logistic
##                                       "t" stands for t distirbution
##         start: start: the starting value for optim, can be NULL
################################################################
M_estimate = function(X, family, m, start = NULL){
  
  gfun = generate_gfun(family)
  
  R = apply(X,2,rank) ## compute rank for the bivariate data
  
  optimal = function(par){optim_fun(par,family,gfun,m,Rank = R)}
  
  if(family == "alog"){
    
    if(is.null(start)) start = c(cor(X[,1], X[,2]), 0.5, 0.5) 
    
    mini = optim(par = start, optimal, 
                  method = "L-BFGS-B", lower = c(0,0,0), upper = c(1,1,1), 
                  control = list(maxit = 1000, factr = 1e10))}
  
  if(family == "t"){
    
    if(is.null(start)) start = c(2,cor(X[,1], X[,2])) 
    
    mini = optim(par = start,optimal,
                  method = "L-BFGS-B", lower = c(1,0.05), upper = c(10,0.99), 
                  control = list(maxit = 1000, factr = 1e10))}
  
  if(family == "hr"){
    
    mini = optim(0, optimal, method = "Brent", lower = 0, upper = 10,
                  control = list(maxit = 10000, abstol = 1e-20, reltol = 1e-20))}
  
  
  if(family == "log"){
    
    mini = optim(0, optimal, method = "Brent", lower = 0, upper = 1)}
  
  return(mini$par)
}

################################################################
#### Hill estimator for tail index
## This function estimate the tail index
## Inputs: dat: a data vector
##         k: sample fraction for hill estimator
################################################################
tail_estimate = function(dat, k){
  n = length(dat)
  if(k> length(which(dat>0))) return(c(k, NA))
  x = sort(dat)
  y = x[(n-k+1):n]
  tk = sum(log(y)-log(x[n-k]))/k
  return(tk)
}

################################################################
#### Nonparametric estimation of quantile Q_Y
## This function estimate VaR nonparametrically
## Inputs: dat: a data vector
##         k: sample fraction for quantile estimation
##         p: the level of quantile
##         gamma: tail index gamma
################################################################
Qvar_esti = function(dat,gamma,p,k){
  n = length(dat)
  X = sort(dat)[n-k]
  Qvar=X*(k/(n*p))^(gamma)
  return(Qvar)
}

################################################################
#### eta_star estimation
## This function estimate eta_star by solving R(1,eta) = p
## Inputs: par_hat: estimates of paramters in TD function
##         p: a bivariate vector of the risk level
##         (the first for VaR, the second for CoVaR)
##         family: name of distribution
################################################################
eta_estimate = function(par_hat, p, family){
  c = p[2]/p[1]
  fun_0 = function(eta){ 
    tail_dependence(X = c(1,eta), par = par_hat, family) - p[2]}
  
  eta_hat = uniroot(fun_0,interval = c(0.0001,c))$root
  return(eta_hat/c)
}

################################################################
#### Main code for CoVaR estimation
## This function estimate CoVaR with the proposed method
## Inputs: Data: two-dimensional dataframe
##         (the first column is institution (conditional event), the second is market)
##         par_hat: estimators of parameters in TD function;
##         if NULL, the function will do the M-estimation befor estimate CoVaR
##         group: a single or a vector of distributions that considered for M-estimation
##         m: a vector of the sample fractions for M-estimator corresponding to group (if par_hat is Null)
##         k: a vector of sample fraction for CoVaR estimation of market index
##         (the first for tail index, the second for VaR estimation)
##         p: a bivariate vector of the risk level
##         (the first for VaR, the second for CoVaR)
################################################################
CoVaR_est = function(Data, par_hat = NULL, group = NULL, m = NULL, k, p){
  
  eta_est = rep(0, length(group))  ### a vector of eat estimates
  
  M_est = list()  ### a list to record M-estimates for different family groups
  
  if(!is.null(par_hat)) M_est = par_hat
  
  if(is.null(par_hat)){
    
    if(length(m)!= length(group)) stop("`m` and `group` must be the same length")
    
    for(i in 1:length(group)){
      
      M_est[[i]] = M_estimate(Data, family = group[i], m = m[i])
      
    }}
  
  for(i in 1:length(group)){
    
    eta_est[i] = eta_estimate(par_hat = M_est[[i]], p = p, family = group[i])
    
  }
  
  sys_tail = tail_estimate(dat = Data[,2], k=k[1])
  
  VaR_est = Qvar_esti(dat = Data[,2], gamma = sys_tail, 
                      k = k[2], p = p[2])
  
  covar_hat = VaR_est*(eta_est)^(-sys_tail)
  
  return(covar_hat)
}



################################################################



################################################################
#### data filtering
## This function filter the original time series data
## Inputs: dat: a data vector
##         (the first column is institution (conditional event), the second is market)
##         forecast: a logical, indicting one-step ahead forecasting or not
##         n_out: the length of one-step ahead forecasting, only value when forecast = TRUE
################################################################
filtering = function(dat, forecast = TRUE, n_out = 1, model = "sGARCH"){
  
  spec = ugarchspec(variance.model = list(model = model, garchOrder = c(1,1)), 
                     mean.model=list(armaOrder=c(1,0),include.mean=T),
                     distribution.model="sstd")
  
  if(forecast){
    
    mean_forecast = sigma_forecast = rep(0, n_out)
    
    for(i in 1:n_out){
      fit = ugarchfit(spec = spec, data = dat[i:(length(dat) - n_out - 1 + i)], 
                       solver.control=list(trace = 1),solver = "hybrid")
      ### rolling window forecasting, every time use length(dat) - n_out samples to fit.
      
      forc = ugarchforecast(fit, n.ahead = 1)
      
      slotNames(fit)
      slotNames(forc)
      mean_forecast[i] = as.numeric(forc@forecast$seriesFor)
      sigma_forecast[i] = as.numeric(forc@forecast$sigmaFor)
    }
    
    fit = ugarchfit(spec = spec, data = dat[1:(length(dat) - n_out)],
                     solver.control=list(trace = 1),solver = "hybrid")
    resid =  array(residuals(fit))
    sigma = array(sigma(fit))
    coef_est = coef(fit)
    
    return(list(residuals = resid/sigma, coef = coef_est,
                mean.forecast = mean_forecast, sigma.forecast = sigma_forecast))}
  
  if(!forecast){
    
    fit = ugarchfit(spec = spec, data = dat,solver.control=list(trace = 1),solver = "hybrid")
    resid =  array(residuals(fit))
    sigma = array(sigma(fit))
    return(resid/sigma)
  }
}


################################################################
#### score function
## This function calculate score
## Inputs: x: verifying observation
##         r: forecast
##         level: upper risk level
################################################################
score_fun = function(x, r, level){
  s = (level-ifelse(x>r, 1,0))*r + ifelse(x>r, 1,0)*x
  return(s)
}

################################################################
#### Main code to calculate average score
## This function calculate average scores
## Inputs: forecast_CoVaR: a (n*d) matrix of CoVaR forecasts for different parametric models
##                         d is the number of models; n is the length of forecasts
##         forecast_VaR: a vector of VaR forecasts for institution
##         test: a two-dimentional dataframe (institution, market)
##         level: upper risk level for CoVaR
################################################################
average_score = function(forecast_CoVaR, forecast_VaR, test, level){
  ## the length of forecast_VaR and test should be the same as nrow(forecast_VaR)
  ind = which(test[,1]>forecast_VaR)
  subsample = test[ind,2]
  forecast_CoVaR = forecast_CoVaR[ind,]
  ave_score = rep(0, ncol(forecast_CoVaR))
  
  for(j in 1:ncol(forecast_CoVaR)){
   est_CoVaR = forecast_CoVaR[,j]
   score = score_fun(x = subsample, r = est_CoVaR,level = level)
   ave_score[j] = mean(score)
  }
  return(ave_score)
}

################################################################
#### unconditional test
## This function conduct unconditional test
## Inputs: forecast_CoVaR: a (n*d) matrix of CoVaR forecasts for different parametric models/methods
##                         d is the number of methods; n is the length of forecasts
##         forecast_VaR: a vector of VaR forecasts for institution
##         test??a two-dimentional dataframe (institution, market)
##         level: risk level for VaR and CoVaR
################################################################
unconditional_test = function(forecast_CoVaR, forecast_VaR, test, level){
  n = nrow(test)
  indi_j = ifelse(test[,1] > forecast_VaR, 1, 0)
  n_VaR = sum(indi_j)
  given.p = level[2]
  ind = which(test[,1] > forecast_VaR)
  test_new = test[ind,2]
  forecast_CoVaR = forecast_CoVaR[ind,]
  re = matrix(0, ncol = 2, nrow = ncol(forecast_CoVaR))
  for(j in 1:ncol(forecast_CoVaR)){
  indi_js = ifelse(test_new > forecast_CoVaR[,j], 1, 0)
  n1_CoVaR = sum(indi_js) ## number of violations
  if(n1_CoVaR == 0){
      logL0 = n1_CoVaR*log(given.p) + (length(test_new)-n1_CoVaR)* log(1-given.p) ## likelihood under H0
      logL1 = n1_CoVaR*log(n1_CoVaR/length(test_new)) + (length(test_new)-n1_CoVaR)*log(1-n1_CoVaR/length(test_new)) ## likelihood under H1
      test_stat = 2*(logL1-logL0)
    }
  if(n1_CoVaR!=0){
      logL0 = n1_CoVaR*log(given.p) + (length(test_new)-n1_CoVaR)* log(1-given.p) ## likelihood under H0
      logL1 = n1_CoVaR*log(n1_CoVaR/length(test_new)) + (length(test_new)-n1_CoVaR)*log(1-n1_CoVaR/length(test_new)) ## likelihood under H1
      test_stat = 2*(logL1-logL0)}
  p_value = pchisq(q = test_stat, df = 1, lower.tail = FALSE)
  re[j,] = c(n1_CoVaR, p_value)
  }
  test_stat_VaR = 2*(n_VaR*log(n_VaR/n) + (n-n_VaR)* log(1-n_VaR/n)) - 
    2*(n_VaR*log(level[1]) + (n-n_VaR)* log(1-level[1]))
  pvalue_VaR = pchisq(q = test_stat_VaR, df = 1, lower.tail = FALSE)
  re = rbind(c(n_VaR, pvalue_VaR), re)
  colnames(re) = c("exceedance", "p_value")
  return(re)
}


################################################################
## Code for other two method

################################################################
#### CoVaR estimation with method in Nolde and Zhang
## This function is used to calculate CoVaR with method in Nolde and Zhang
## Input: fit is the return from New_cov_update
##        q is quantile
##        VaR is VaR for individual stock 
################################################################
CoVaR_EVT = function(x, q, fit, VaR){ # x is CoVaR
  z = c(VaR, x)
  diag_inv_EVT = matrix(c(1/fit$EVT[5], 0,0,1/fit$EVT[6]),2,2)
  tmp_value_EVT = diag_inv_EVT %*% z
  z_EVT =  tmp_value_EVT[2] / tmp_value_EVT[1] - fit$EVT[3]
  theorm_EVT = 1-Thrm4_1(z_EVT, fit$EVT[1], fit$EVT[2], fit$EVT[3], fit$EVT[4])
  return( (theorm_EVT - q)^2*10^10 )
}

################################################################
#### Calculate CoVaR for Skew-T method
## This function is used to calculate CoVaR with skew-t method in Girardi2013
## Input: VaR is VaR for individual stock 
##        rho: correlation parameter in skew-t distribution
##        alpha: skew parameter in skew-t distribution
##        nu: degree of freedom
##        p: a vector of risk level
################################################################
CoVaR_Skew = function (x, VaR, xi, Omega, alpha, nu, p) {
  joint = 1- sn::pmst(c(VaR, Inf), xi = xi, Omega = Omega, alpha = alpha, nu = nu) -
    sn::pmst(c(Inf, x), xi = xi, Omega = Omega, alpha = alpha, nu = nu) + 
    sn::pmst(c(VaR, x), xi = xi, Omega = Omega, alpha = alpha, nu = nu)
  return( (joint - p[1]*p[2])^2*10^10 )
}

################################################################
#### Main code to estimate CoVaR with other two method
## This function is used to calculate CoVaR with skew-t method in Girardi2013
## Input: Data: two-dimensional dataframe
##        fit_par: fitted coefficents from ugarchfit of institution
##        VaR: VaR forecast of institution
##        k: if VaR = NULL, k is the sample fraction vector to estimate VaR of institution
##        level: a vector of risk levels
################################################################
CoVaR_FP = function(Data, fit_par = NULL, VaR = NULL, level){
  if(is.null(VaR)){
    nu=fit_par["shape"]; ga=fit_par["skew"]
    VaR = qdist("sstd",p=1-level[1],mu=0,sigma=1,shape=nu,skew=ga)
  }
  
  ##### Parametric method ###########
  fit_skew = mst.mle(y = Data)$dp
  
  tmp_Skew=as.numeric(nlminb(start=4, objective=CoVaR_Skew,VaR = VaR, 
                                  xi = fit_skew$beta, Omega = fit_skew$Omega,
                                  alpha = fit_skew$alpha, nu = fit_skew$df,
                                  p = level)$par)
  
  return(tmp_Skew)
}

CoVaR_NZ = function(Data, VaR = NULL, k = NULL, level){
  tmp = NA
  
  if(is.null(VaR)){
    VaR = VaR_estimate(dat = Data[,1], k, p = level[1])
  }
  
  fit =  New_cov_update(Data)
  
  tmp = as.numeric(nlminb(start=4, objective=CoVaR_EVT,q = level[2], 
                               fit = fit, VaR = as.numeric(VaR))$par)
  
  return(tmp)
}




################################################################
## Below are functions used in EVT method of Nolde and Zhang
## Copied from Nolde and Zhang


################################################################
##### Functions used to calculate theoretical value in Theorem 3.2
## alpha1 and alpha2 are shape parameters
## rho is off-diagonal element in standardized scale matrix
## v (or nu) is tail index 
## xi is location parameter
############################################################

############################################################
##### K(z;alpha1, alpha2, rho ,v)
############################################################
#### Upper part of lambda function: T_1 * Q
int_up = function(x, alpha1, alpha2, rho, v){
  x1 = x[1]
  x2 = x[2]
  Q_up = ((x1)^2 - 2*rho*(x1)*(x2) + (x2)^2 )/(1-rho^2)
  T_up = pt( (alpha1*x1 + alpha2*x2) * sqrt((v+2)/ Q_up), v+2)
  Q_up^(-(v+2)/2)*T_up
}
#### Lower part of lambda function: Q
int_down = function(x, alpha1, alpha2, rho, v){
  x1 = x[1]
  x2 = x[2]
  Q_up = ((x1)^2 - 2*rho*(x1)*(x2) + (x2)^2 )/(1-rho^2)
  Q_up^(-(v+2)/2)
}
#### Numerical function of K(z, alpha1, alpha2, rho, v)
K = function(z, alpha1, alpha2, rho, v){
  K_up  = 2 * hcubature(f=int_up,lowerLimit = c(1, z), upperLimit = c(10^6,10^6), tol = 1e-8,  alpha1 = alpha1, alpha2 = alpha2, rho = rho, v = v)$integral / hcubature(f=int_down,lowerLimit = c(1, z), upperLimit = c(10^6,10^6), tol = 1e-8,  alpha1 = alpha1, alpha2 = alpha2, rho = rho, v = v)$integral
  alpha1_bar = (alpha1 + rho *alpha2) / sqrt(1 + alpha2^2 * (1-rho^2))
  K_down = 2 * pt(alpha1_bar * sqrt(v+1), df = v+1)
  return(K_up/K_down)
}

##############################################################
#### Calculate the theoretic value of Pr(Y<y|X>x) assuming (X,Y) are bivariate skew-t distribution. 
##############################################################
Theory_value = function(x, xi, alpha, rho, v){
  Omega = matrix(c(1,rho,rho,1),2,2)
  Pxy = pmst(c(Inf,x[2]), xi , Omega, alpha, v, abseps = 1e-14, maxpts = 1000000)  - pmst(c(x[1],x[2]), xi , Omega, alpha, v,abseps = 1e-14, maxpts = 100000)
  Px =  1 - pmst(c(x[1],Inf), xi , Omega, alpha, v,abseps = 1e-14, maxpts = 1000000)
  Pxy/Px
}
Theory_value_omega = function(x, xi, alpha, Omega, v){
  Pxy = pmst(c(Inf,x[2]), xi , Omega, alpha, v)  - pmst(c(x[1],x[2]), xi , Omega, alpha, v)
  Px =  1 - pmst(c(x[1],Inf), xi , Omega, alpha, v)
  Pxy/Px
}

##############################################################
#### Calculate the empirical value of Pr(Y<y|X>x)
##############################################################
Empirical_value = function(data, quan_x,quan_y){
  n = nrow(data)
  tmp = subset(data, data[,1]> quan_x)
  Px = nrow(tmp)/n
  Pxy = nrow(subset(tmp, tmp[,2] > quan_y))/n
  1 - Pxy/Px
}

##############################################################
#### Calculate the simulated value of Pr(Y<y|X>x) when y = rho x
#### based on the Theorem 3.2 in the paper (it was Theorem 4.1)
################################################################
Thrm4_1 = function(z, alpha1, alpha2,rho, v){
  k = K( rho+z, alpha1, alpha2, rho, v)
  sign = sign(rho+z)
  T1_bar1 = 1 - pt(z * sqrt(v+1)/sqrt(1-rho^2), df = v+1)
  T1_bar2 = 1 - pt(sign * sqrt(v+1)/sqrt(1-rho^2) * (1/(rho+z) - rho), df = v+1)
  1 - k * (T1_bar1 + sign/(abs(rho+z))^v * T1_bar2 )
}

##############################################################
#### psi function in equation (4.1) defined in paper <Estimation of extreme risk regions #### under multivariate regular variation>
##############################################################
func = function(x,rho,v){ #Function A(theta)
  y = (1 + rho^2 * cos(2*x) + rho * sqrt(1-rho^2)*sin(2*x) )^(v/2)
  return(y)
}
func_cov = function(x,rho,sd1, sd2,v){
  y = ( sd1^2 * cos(x)^2 + sd2^2 * (sin(x)^2 + rho^2 * cos(2*x) + rho * sqrt(1-rho^2)*sin(2*x) ) )^(v/2)
  return(y)
}
q_x = function(x1,x2, alpha1, alpha2, rho, v){
  x = c(x1,x2)
  Omega = matrix(c(1,rho,rho,1),2,2)
  Q_up = ((x1)^2 - 2*rho*(x1)*(x2) + (x2)^2 ) /(1-rho^2)
  T_up = pt( (alpha1*x1 + alpha2*x2) * sqrt( (v+2)/ Q_up), v+2)
  tmp = det(Omega)^(-0.5)
  #as.numeric(Q_up^(-(v+2)/2) * T_up * tmp / pi * v)
  int = integrate(func,0,2*pi, rho = rho, v= v)$value
  as.numeric(Q_up^(-(v+2)/2) * T_up * tmp *2 * v / int)
}

q_x_cov = function(x1,x2, alpha1, alpha2, rho, v, sd1, sd2){
  x = c(x1,x2)
  Omega = matrix(c(sd1^2,rho*sd1 * sd2,rho*sd1 * sd2,sd2^2),2,2)
  Q_up = (sd2^2 * (x1)^2 - 2*rho*(x1)*(x2) * sd1 * sd2 + sd1^2 * (x2)^2 ) /( sd1^2 * sd2^2*(1-rho^2) )
  T_up = pt( (alpha1*x1/sd1 + alpha2*x2/sd2) * sqrt( (v+2)/ Q_up), v+2)
  tmp = det(Omega)^(-0.5)
  #as.numeric(Q_up^(-(v+2)/2) * T_up * tmp / pi * v)
  int = integrate(func_cov,0,2*pi, rho = rho, sd1 = sd1, sd2 = sd2, v= v)$value
  as.numeric(Q_up^(-(v+2)/2) * T_up * tmp *2 * v / int)
}


###########################################################
######## Starting from here, the functions are used to do inference as described in Section 4 
## alpha1 and alpha2 are shape parameters
## rho is off-diagonal element in standardized scale matrix
## v (or nu) is tail index 
## xi is location parameter
############################################################

##############################################################
#### spectral measure functions to draw graph between 0 and 2*pi
##############################################################
q_theta = function(x){
  x1 = cos(x)
  x2 = sin(x)
  x = c(x1,x2)
  Q_up = ((x1)^2 - 2*rho*(x1)*(x2) + (x2)^2 ) /(1-rho^2)
  T_up = pt( (alpha1*x1 + alpha2*x2) * sqrt( (v+2)/ Q_up), v+2)
  tmp = det(Omega)^(-0.5)
  # as.numeric(Q_up^(-(v+2)/2) * T_up * tmp / pi)
  int = integrate(func,0,2*pi, rho = rho, v= v)$value
  as.numeric(Q_up^(-(v+2)/2) * T_up * tmp *2 / int)
}


##############################################################
#### Negative log likelihood function (4.1) to estimate alpha1 alph2 and rho
##############################################################
neglogfn = function(param, sim_y, v, prob = 0.80){
  alpha1 = param[1]
  alpha2 = param[2]
  rho = param[3]
  if(rho < -1 | rho > 1){ return(1.e10)}
  if(abs(alpha1)>20){return(1.e10)}
  if(abs(alpha2)>20){return(1.e10)}
  sim_y_iter = sim_y
  n = nrow(sim_y)
  p = c(rep(NA,n))
  w = sim_y
  for(i in 1:n){
    p[i] = sqrt(sim_y_iter[i,1]^2 + sim_y_iter[i,2]^2)
    w[i,] = sim_y_iter[i,]/p[i]
  }
  p = cbind(p,w)
  colnames(p) = c("r", "cos", "sin")
  p = data.frame(p)
  
  sort.p = p[order(p$r), ]
  cut = prob* n
  cut_point = sort.p[cut:(n-2),]
  
  value = -sum(log(1/v*q_x(cut_point[,2],cut_point[,3], alpha1, alpha2, rho, v)))
  return(value)
}

neglogfn1 = function(param, sim_y,v, prob = 0.80){
  alpha1 = param[1]
  alpha2 = param[2]
  rho = param[3]
  sd1 = param[4]
  sd2 = param[5]
  if(rho < -1 | rho > 1){ return(1.e10)}
  if(abs(alpha1)>20){return(1.e10)}
  if(abs(alpha2)>20){return(1.e10)}
  if(sd1<0){return(1.e10)}
  if(sd2<0){return(1.e10)}
  n = nrow(sim_y)
  n = nrow(sim_y)
  p = c(rep(NA,n))
  w = sim_y
  for(i in 1:n){
    p[i] = sqrt(sim_y[i,1]^2 + sim_y[i,2]^2)
    w[i,] = sim_y[i,]/p[i]
  }
  p = cbind(p,w)
  colnames(p) = c("r", "cos", "sin")
  p = data.frame(p)
  
  sort.p = p[order(p$r), ]
  cut = prob * n
  cut_point = sort.p[cut:(n),]
  
  value = -sum(log(1/v*q_x_cov(cut_point[,2],cut_point[,3], alpha1, alpha2, rho, v, sd1, sd2)))
  return(value)
}


###################################################
#-----------> function Hill to estimate v
# paper: A Multivariate Hill Estimator
###################################################
Hill_small = function(sim_y, start, end){
  ncol = ncol(sim_y)
  n = nrow(sim_y)
  if(ncol > 0){
    p = c(rep(NA,n))
    for(i in 1:n){
      p[i] = sqrt(sum(sim_y[i,]^2))
    }
    v_temp = hill(p, start = start,end = end, option = "xi")
  }
  #if(ncol == 1){
  #  v_temp = hill(sim_y, start = start,end = end, option = "xi")
  #}
  gamma = v_temp$y
  Num = end-start +1
  
  W = Diagonal(n=Num, x=sqrt(v_temp$x))
  Z = cbind((rep(1,Num)),v_temp$x)
  b = solve(t(Z)%*%t(W)%*%W%*%Z) %*% t(Z)%*%t(W)%*%W%*%gamma
  est_v = 1/b[1]
  return(est_v)
}

#### Plot f_theta
f_theta_plot = function(theta){
  L = matrix(c(1,0,rho,sqrt(1-rho^2)),2,2)
  alpha_tmp = L %*% c(alpha1,alpha2)
  alpha1_star = alpha_tmp[1]
  alpha2_star = alpha_tmp[2]
  pt(sqrt(2) * (alpha1_star * cos(theta) + alpha2_star * sin(theta)), df=2)/pi
}
#### f_theta function
f_theta = function(rho, theta, alpha1, alpha2){
  L = matrix(c(1,0,rho,sqrt(1-rho^2)),2,2)
  alpha_tmp = L %*% c(alpha1,alpha2)
  alpha1_star = alpha_tmp[1]
  alpha2_star = alpha_tmp[2]
  pt(sqrt(2) * (alpha1_star * cos(theta) + alpha2_star * sin(theta)), df=2)/pi
}




##################################################################
## Main code to do inference described in section 4. 
## New_cov_update is the most general case.
## Input:  sim_y, a matrix with two vectors
## Output: EVT -- estimation after re-estimating center parameter xi
##         EVT_old -- estimation without re-estimating ceter parameter xi, using Algorithm 1
##         xi -- estimation of center parameter
##         alpha -- estimation of skew parameters from Algo 2 (only for comparision)
##         sd1, sd2, rho -- estimation of scale parameter from Algo 2 (only for comparision)
##         v -- estimation of the tail index
##################################################################
New_cov_update = function(sim_y, method="combine", prob = 0.85){
  
  ## Initial value
  set.seed(123)
  CovRob = cov.rob(sim_y, cor = TRUE)
  
  xi_est = CovRob$center
  est_rho = CovRob$cor[1,2]
  cov_est = CovRob$cov
  
  n = nrow(sim_y)
  mxi_up = as.matrix(xi_est)
  mxi_up = t(mxi_up[,rep(1,n)])
  sim_y_iter = sim_y - mxi_up
  
  #### Extreme value theory
  est_v = Hill_small(sim_y_iter, 5, as.integer(0.5*n))
  fit = nlm(neglogfn1, c(-1,1, est_rho, sqrt(cov_est[1,1]),sqrt(cov_est[2,2])),  sim_y = sim_y_iter, v = est_v, prob = prob)
  while (!(fit$code %in% c(1,2))){
    fit = nlm(neglogfn1, c(-runif(1), runif(1), est_rho, sqrt(cov_est[1,1]),sqrt(cov_est[2,2])),  sim_y = sim_y_iter, v = est_v, prob = prob)
  }
  
  est_alpha1  =  fit$estimate[1]
  est_alpha2 = fit$estimate[2]
  est_rho = fit$estimate[3]
  est_sd1 = fit$estimate[4]
  est_sd2 = fit$estimate[5]
  est_EVT_old = c(est_alpha1,est_alpha2,est_rho,est_v, est_sd1, est_sd2)
  
  #####################
  ### Algorithm 3######
  #####################
  fit_st1 = mst.fit(y = sim_y, plot.it=F)
  xi_est = as.numeric(fit_st1$dp$beta)
  cov_est = fit_st1$dp$Omega
  rho_est = cov_est[1,2]/sqrt(cov_est[1,1]*cov_est[2,2])
  alpha1_est = fit_st1$dp$alpha[1]
  alpha2_est = fit_st1$dp$alpha[2]
  
  ### Only adopt estimated xi from above algorithm
  ### re-update EVT
  mxi_up = as.matrix(xi_est)
  mxi_up = t(mxi_up[,rep(1,n)])
  sim_y_iter = sim_y - mxi_up
  
  #### Extreme value theory
  fit = nlm(neglogfn1, c(alpha1_est, alpha2_est, rho_est, sqrt(cov_est[1,1]),sqrt(cov_est[2,2])),  sim_y = sim_y_iter, v = est_v, prob = prob)
  while (!(fit$code %in% c(1,2))){
    fit = nlm(neglogfn1, c(-runif(1)+1, runif(1), rho_est, sqrt(cov_est[1,1]),sqrt(cov_est[2,2])),  sim_y = sim_y_iter, v = est_v, prob = prob)
  }
  
  est_alpha1 =  fit$estimate[1]
  est_alpha2 = fit$estimate[2]
  est_rho = fit$estimate[3]
  est_sd1 = fit$estimate[4]
  est_sd2 = fit$estimate[5]
  
  est_EVT = c(est_alpha1,est_alpha2,est_rho,est_v, est_sd1, est_sd2)
  return(list(est_v = fit_st1$dp$df, EVT = est_EVT, EVT_old = est_EVT_old))
}

##################################################################
## Main code to do inference described in section 4. 
## New_update is for the case that the scale matrix is standard (no need to estimate parameter sd), and the locatoin parameter is c(0,0)
## Input:  sim_y, a matrix with two vectors
## Output: EVT -- estimation after re-estimating center parameter xi
##         EVT_old -- estimation without re-estimating ceter parameter xi, using Algorithm 1
##         xi -- estimation of center parameter
##         alpha -- estimation of skew parameters from Algo 2 (only for comparision)
##         rho -- estimation of scale parameter from Algo 2 (only for comparision)
##         v -- estimation of the tail index
##################################################################
New_update = function(sim_y, method="combine", prob = 0.85){
  #### Extreme value theory
  set.seed(123)
  CovRob = cov.rob(sim_y, cor = TRUE)
  fit_st1 = mst.fit(y = sim_y, plot.it=F)
  xi_est = as.numeric(fit_st1$dp$beta)
  
  # xi_est = CovRob$center
  rho_est = CovRob$cor[1,2]
  
  n = nrow(sim_y)
  mxi_up = as.matrix(xi_est)
  mxi_up = t(mxi_up[,rep(1,n)])
  sim_y_iter = sim_y - mxi_up
  
  #### Extreme value theory
  est_v = Hill_small(sim_y_iter, as.integer(0.01*n), as.integer(0.25*n))
  fit = nlm(neglogfn, c(-1,0, rho_est),  sim_y = sim_y_iter, v = est_v, prob = prob)
  est_EVT_old = c(fit$estimate[1], fit$estimate[2],fit$estimate[3],est_v)
  
  #####################
  ### Algorithm 3######
  #####################
  cov_est = fit_st1$dp$Omega
  rho_est = cov_est[1,2]/sqrt(cov_est[1,1]*cov_est[2,2])
  alpha1_est = fit_st1$dp$alpha[1]
  alpha2_est = fit_st1$dp$alpha[2]
  
  #### Extreme value theory
  fit1 = nlm(neglogfn, c(alpha1_est,alpha2_est,rho_est),  sim_y = sim_y_iter, v = est_v, prob = prob)
  
  est_EVT = c(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3],est_v)
  return(list(est_v = est_v, EVT = est_EVT, EVT_old = est_EVT_old))
}


##################################################################
##### FUNCTIONS FOR BACKTESTING
##################################################################

##################################################################
##### TLM test 
## The function computes TLM's for three significance levels (1%, 5%, 10%)
## comparison of two methods only
## Inputs: tmp: list of scores for one method ("proposed" method)
##         tmpS: list of scores for one method for the other method (to be compared to, "standard" method)
##         cor: HAC estimator for the asymptotic variance of the average relative scores (if False, put correlation as 0)
##         p: risk level for CoVaR
##################################################################
pool_efp.test = function(s1,s2,type="two-sided", cor = TRUE){
  ### s1, s2 should be a list, each list component stands for a firm; s2 is the reference method
  
  # differences in scores
  p = length(s1)
  # HAC estimator for the asymptotic variance of the average relative scores
  # HAC: heteroskedasticity and autocorrelation-consistent variance estimator
  # using Parzen window as the lag window (set of weights)
  # M: trancation point, set to approx. 2sqrt(n)
  scale_d = list()
  for(i in 1:p){
    d = s1[[i]]-s2[[i]]
    n = length(d)
    
    if(!cor){sn = sqrt(sum(d^2)/n)} # assumes zero correlation
    
    if(cor){
      m = ceiling(2*sqrt(n))
      gam = acf(d,lag.max=m,type="covariance",plot=F,na.action = na.fail)$acf
      k1 = 1:ceiling(m/2)
      k2 = (ceiling(m/2)+1):m
      lam = c(1, 2*(1-6*(k1/m)^2+6*(k1/m)^3),2*2*(1-k2/m)^3)
      sn = sqrt(gam %*% lam)}
    
    scale_d[[i]] = d/as.numeric(sn)
  }
  
  # test statistic
  tn = mean(unlist(scale_d))*sqrt(length(unlist(scale_d)))
  
  if (type == "two-sided")
  {
    pv=2*ifelse(tn<0,pnorm(tn),(1-pnorm(tn))) #p-value for two-sided hypothesis
  }
  if (type == "one-sided-ge")
  {
    pv= pnorm(tn)
  }
  if (type == "one-sided-le")
  {
    pv= 1-pnorm(tn)
  }
  # cat("mean=", mean(d),"; var=", var(d), ";sn=",sn, "\n")
  return(list(tn=tn,pvalue=pv,sn=sn))
}

TLMfn.biv = function(tmp, tmpS, cor = TRUE){
  test = pool_efp.test(tmp,tmpS,type="one-sided-ge", cor = cor) # H0^+
  Tn = test$tn ## value of the test statistic
  pvmat = test$pvalue
  pvmatL = pool_efp.test(tmp,tmpS,type="one-sided-le", cor = cor)$pvalue # H0^-
  
  TLM01 = (pvmat = 0.01) + (pvmatL > 0.01)
  TLM05 = (pvmat = 0.05) + (pvmatL > 0.05)
  TLM10 = (pvmat = 0.10) + (pvmatL > 0.10)
  return(c(TLM01=TLM01, TLM05=TLM05, TLM10=TLM10,Tn=Tn))
}

comparative.test = function(VaR_fr, CoVaR_fr, CoVaR_fr_r, test_ins, test_sys, p){
  ### VaR_fr is the forecasted value of VaR of intitutions
  ### CoVaR_fr is forecasted value of CoVaR 
  ### CoVaR_fr_r is forecasted value of CoVaR of reference method
  d = ncol(test_ins)
  n = nrow(test_ins)
  
  tmp = list()
  tmpS = list()
  
  for(i in 1:d){
    ind_pick = which(test_ins[,i] > VaR_fr[,i])
    test_sub = test_sys[ind_pick]
    CoVaR_fr_sub = CoVaR_fr[ind_pick,i]
    CoVaR_fr_r_sub = CoVaR_fr_r[ind_pick,i]
    tmp[[i]] = score_fun(x = test_sub, r = CoVaR_fr_sub, level = p)
    tmpS[[i]] = score_fun(x = test_sub, r = CoVaR_fr_r_sub, level = p)
  }
  re = TLMfn.biv(tmp, tmpS, cor = TRUE)
  
  return(re)
}

###########################################################
#Skew Functions --
# The left functoins are used to fit multivariate skew-t distribution
# which are copied from old library sn.
############################################################

mst.fit = function(X, y, freq, start, fixed.df=NA, plot.it=TRUE,
                    trace=FALSE, ...)
{
  y.name = deparse(substitute(y))
  y.names= dimnames(y)[[2]]
  y = as.matrix(y)
  d = ncol(y)
  if(is.null(d)) d= 1
  if(d>1){
    if(length(y.names)==0){
      dimnames(y) =
        list(dimnames(y)[[1]], outer("V",as.character(1:d),paste,sep=""))
      y.names= as.vector(dimnames(y)[[2]])
    }}
  else
    colnames(y)=y.name
  if(missing(freq)) freq = rep(1,nrow(y))
  n = sum(freq)
  if(missing(X)) {
    X = rep(1,nrow(y))
    missing.X = TRUE }
  else
    missing.X = FALSE
  X   = as.matrix(X)
  qrX = qr(X)
  m   = ncol(X)
  mle = mst.mle(X=X, y=y, freq=freq,  fixed.df=fixed.df, start=start,
                 trace=trace, ...)
  mle$call = match.call()
  beta  = mle$dp$beta
  Omega = mle$dp$Omega
  alpha = mle$dp$alpha
  omega = sqrt(diag(Omega))
  df    = mle$dp$df
  xi    = X %*% beta
  if(plot.it & all(freq==rep(1,length(y)))) {
    if(missing.X) {
      y0  =y
      xi0 = apply(xi,2,mean)}
    else  {
      y0  = y-xi
      xi0 = rep(0,d)
    }
    if(d>1) {
      opt=options()
      options(warn=-1)
      pairs(y0, labels=y.names,
            panel=function(x, y, Y, y.names, xi, Omega, alpha)
            {
              for(i in 1:length(alpha)){
                if(all(Y[,i]==x)) Ix=i
                if(all(Y[,i]==y)) Iy=i
              }
              points(x,y)
              marg = msn.marginal(xi, Omega ,alpha, c(Ix,Iy))
              xi.marg = marg$xi
              Omega.marg = marg$Omega
              alpha.marg = marg$alpha
              x1 = seq(min(x), max(x), length=30)
              x2 = seq(min(y), max(y), length=30)
              dst2.plot(x1, x2, xi.marg, Omega.marg, alpha.marg, df,
                        add=TRUE, col=2)
            },  # end "panel" function
            Y=y0, y.names=y.names, xi=xi0, Omega=Omega, alpha=alpha)
      options(opt)
    }
    else{ # plot for case d=1
      y0=as.vector(y0)
      x=seq(min(pretty(y0,10)),max(pretty(y0,10)),length=100)
      if(missing.X){
        dp0=c(xi0,omega,alpha,df)
        xlab=y.name}
      else {
        dp0=c(0,omega,alpha,df)
        xlab = "residuals"}
      hist(y0, prob=TRUE,  breaks="FD", xlab=xlab, ylab="density", main="")
      lines(x, dst(x,dp0[1],dp0[2],dp0[3],dp0[4]),  col=2)
      if(length(y)<101) points(y0, rep(0,n), pch=1)
      title(y.name)
    }
    cat("Press <Enter> to continue..."); readline()
    par(mfrow=c(1,2))
    pp  = d * qf((1:n)/(n+1),d,df)
    pp2 = qchisq((1:n)/(n+1),d)
    # Xb  = qr.fitted(qrX,y)
    res = qr.resid(qrX,y)
    rad.n  = apply(res    * (res %*% solvePD(var(res))), 1, sum)
    rad.st = apply((y-xi) * ((y-xi) %*% solvePD(Omega)), 1, sum)
    plot(pp2, sort(rad.n), pch=1, ylim=c(0,max(rad.n,rad.st)),
         xlab="Percentiles of chi-square distribution",
         ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for normal distribution", sub=y.name)
    plot(pp, sort(rad.st), pch=1, ylim=c(0,max(rad.n,rad.st)),
         xlab="Percentiles of scaled F distribution",
         ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for skew-t distribution", sub=y.name)
    prob = pf(rad.st/d,d,df)
    mle$mahalanobis = list(distance=rad.st, prob=prob, df=c(d,df))
    cat("Press <Enter> to continue, 'q' to quit...")
    m = readline()
    if(tolower(m) != "q") {
      plot((1:n)/(n+1), sort(pchisq(rad.n,d)), xlab="", ylab="")
      abline(0,1,lty=3)
      title(main="PP-plot for normal distribution", sub=y.name)
      plot((1:n)/(n+1), sort(prob), xlab="", ylab="")
      abline(0,1,lty=3)
      title(main="PP-plot for skew-t distribution", sub=y.name)
    }
    par(mfrow=c(1,1))
    
  } # end ploting
  dev.norm = msn.dev(c(qr.coef(qrX,y),rep(0,d)), as.matrix(X), y, freq)
  test = dev.norm + 2*mle$logL
  p.value =  1-pchisq(test,d+1)
  if(trace) {
    cat("LRT for normality (test-function, p-value): ")
    print(c(test,p.value))
  }
  mle$test.normality = list(LRT=test, df=d+1, p.value=p.value,
                             normal.logL=dev.norm/(-2))
  invisible(mle)
}

#

st.mle = function(X, y, freq,  start, fixed.df=NA, trace=FALSE,
                   algorithm = c("nlminb","Nelder-Mead", "BFGS", "CG", "SANN"),
                   control=list())
{
  y.name  = deparse(substitute(y))
  y = data.matrix(y)
  if(missing(X)) X= matrix(1, nrow=length(y), ncol=1)
  dimnames(y)[[2]] = list(y.name)
  if(missing(start)){
    cp0 = sn.mle(X=X, y=y, plot.it=FALSE, trace=trace)$cp
    m = length(cp0)-2
    cp0[m+2] = cp0[m+2]*0.9
    mle0 = cp.to.dp(cp0)
    start = list(beta=mle0[1:m], Omega=matrix(mle0[m+1]^2,1,1),
                  alpha=mle0[m+2], df=10)
  }
  else {
    m = length(start)-3
    if(m<1) stop("bad start vector")
    start=  list(beta=start[1:m], Omega=matrix(start[m+1]^2,1,1),
                  alpha=start[m+2], df=start[m+3])
  }
  fit = mst.mle(X, y, freq, start=start, fixed.df=fixed.df, trace=trace,
                 algorithm=algorithm, control=control)
  mle = list()
  mle$call= match.call()
  dp = fit$dp
  se = fit$se
  p  = length(dp$beta)
  dp.names = c(if(p==1) "location" else dimnames(dp$beta)[[1]],
                "scale","shape","df")
  mle$dp  = c(dp$beta, sqrt(as.vector(dp$Omega)), dp$alpha, dp$df)
  names(mle$dp) = dp.names
  mle$se = if(all(is.na(se))) NA else
    c(se$beta, mle$dp[p + 1] * se$internal[p + 1],
      se$alpha, dp$df * se$internal[p + 3])
  mle$logL = fit$logL
  mle$algorithm = fit$algorithm
  mle
}


mst.mle = function (X, y, freq, start, fixed.df = NA, trace = FALSE,
                     algorithm = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"),
                     control = list())
{
  algorithm = match.arg(algorithm)
  y.name = deparse(substitute(y))
  y.names = dimnames(y)[[2]]
  y = data.matrix(y)
  X = if (missing(X)) matrix(rep(1, nrow(y)), ncol = 1)
  else data.matrix(X)
  if (missing(freq)) freq = rep(1, nrow(y))
  x.names = dimnames(X)[[2]]
  d = ncol(y)
  n = sum(freq)
  m = ncol(X)
  if (missing(start)) {
    qrX = qr(X)
    beta = as.matrix(qr.coef(qrX, y))
    Omega = matrix(var(qr.resid(qrX, y)), d, d)
    omega = sqrt(diag(Omega))
    alpha = rep(0, d)
    df = ifelse(is.na(fixed.df), 10, fixed.df)
    if (trace) {
      cat("mst.mle: dp=", "\n")
      print(c(beta, Omega, alpha))
      cat("df:", df, "\n")
    }
  }
  else {
    if (!is.na(fixed.df))
      start$df = fixed.df
    if (all(names(start) == c("beta", "Omega", "alpha", "df"))) {
      beta = start$beta
      Omega = start$Omega
      alpha = start$alpha
      df = start$df
    }
    else stop("start parameter is not in the form that I expected")
  }
  eta = alpha/sqrt(diag(Omega))
  Oinv = solvePD(Omega)
  upper = chol(Oinv)
  D = diag(upper)
  A = upper/D
  D = D^2
  if (d > 1)
    param = c(beta, -log(D)/2, A[!lower.tri(A, diag = TRUE)], eta)
  else
    param = c(beta, -log(D)/2, eta)
  if (is.na(fixed.df))
    param = c(param, log(df))
  if(algorithm == "nlminb"){
    opt = nlminb(param, objective = mst.dev, gradient = mst.dev.grad,
                  control = control,  X = X, y = y, freq = freq,
                  trace = trace, fixed.df = fixed.df)
    info = num.deriv2(opt$par, FUN="mst.dev.grad", X=X, y=y,
                       freq=freq, fixed.df = fixed.df)/2
    opt$value =  opt$objective
  }
  else{
    opt = optim(param, fn = mst.dev, gr = mst.dev.grad,
                 method = algorithm, control = control, hessian = TRUE,
                 X = X, y = y, freq = freq, trace = trace, fixed.df = fixed.df)
    info = opt$hessian/2
  }
  dev   = opt$value
  param = opt$par
  opt$name = algorithm
  if (trace) {
    cat("Message from optimization routine:", opt$message, "\n")
    cat("deviance:", dev, "\n")
  }
  beta = matrix(param[1:(m * d)], m, d)
  D = exp(-2 * param[(m * d + 1):(m * d + d)])
  A = diag(d)
  i0 = m*d+d*(d+1)/2
  if(d>1)  A[!lower.tri(A,diag=TRUE)] = param[(m*d+d+1):i0]
  eta = param[(i0 + 1):(i0 + d)]
  if (is.na(fixed.df))
    df = exp(param[i0 + d + 1])
  else df = fixed.df
  Oinv = t(A) %*% diag(D,d,d) %*% A
  Omega = solvePD(Oinv)
  omega = sqrt(diag(Omega))
  alpha = eta * omega
  dimnames(beta) = list(x.names, y.names)
  dimnames(Omega) = list(y.names, y.names)
  if (length(y.names) > 0) names(alpha) = y.names
  if (all(is.finite(info))) {
    qr.info = qr(info)
    info.ok = (qr.info$rank == length(param))
  }
  else info.ok = FALSE
  if (info.ok) {
    se2 = diag(solve(qr.info))
    if (min(se2) < 0)
      se = NA
    else {
      se = sqrt(se2)
      se.beta = matrix(se[1:(m * d)], m, d)
      se.alpha = se[(i0 + 1):(i0 + d)] * omega
      dimnames(se.beta)[2] = list(y.names)
      dimnames(se.beta)[1] = list(x.names)
      names(se.alpha) = y.names
      se.df = df * se[i0 + d + 1]
      se = list(beta = se.beta, alpha = se.alpha, df = se.df,
                 internal = se, info = info)
    }
  }
  else se = NA
  dp = list(beta = beta, Omega = Omega, alpha = alpha, df = df)
  list(call = match.call(), logL = -dev/2, deviance = dev,
       dp = dp, se = se, algorithm = opt)
}


solvePD = function(x){ 
  # inverse of a symmetric positive definite matrix
  u = chol(x, pivot = FALSE)
  if(prod(diag(u)) == 0) 
    stop("matrix not positive definite")
  # ui = backsolve(u,diag(ncol(x)))
  # ui %*% t(ui)
  chol2inv(u)
}


mst.dev = function(param, X, y, freq=rep(1,nrow(X)), fixed.df=NA, trace=FALSE)
{
  # Diag = function(x) diag(x, nrow=length(x), ncol=length(x))
  d = ncol(y)
  # if(missing(freq)) freq=rep(1,nrow(y))
  n = sum(freq)
  m = ncol(X)
  beta=matrix(param[1:(m*d)],m,d)
  D = exp(-2*param[(m*d+1):(m*d+d)])
  i0 = m*d+d*(d+1)/2
  A = diag(d)
  if(d>1) A[!lower.tri(A,diag=TRUE)] = param[(m*d+d+1):i0]
  eta = param[(i0+1):(i0+d)]
  if(is.na(fixed.df)) df = exp(param[i0+d+1])
  else df = fixed.df
  Oinv = t(A) %*% diag(D,d,d) %*% A
  # Omega = solvePD(Oinv)
  u =  y - X %*% beta
  Q = apply((u %*% Oinv)*u,1,sum)
  L = as.vector(u %*% eta)
  logDet= sum(-log(D))
  if(df < 10000)  {
    const= lgamma((df + d)/2)- lgamma(df/2)-0.5*d*logb(df)
    DQ =  (df+d) * sum(freq *logb(1+Q/df))
    L. = L*sqrt((df+d)/(Q+df))
  }
  else {
    const = (-0.5*d*logb(2)+ log1p((d/2)*(d/2-1)/df))
    DQ = if(df<Inf) (df+d) * sum(freq *log1p(Q/df)) else sum(freq*Q)
    L. = L*sqrt((1+d/df)/(1+Q/df))
  }
  dev = (n*(logDet - 2*const+ d*logb(pi)) + DQ
          -2*sum(freq * (log(2)+log.pt(L., df+d))))
  if(trace) cat("mst.dev: ",dev, "\n")
  dev
}


mst.dev.grad = function(param, X, y, freq=rep(1,nrow(X)), fixed.df=NA,
                         trace=FALSE)
{
  # Diag = function(x) diag(x, nrow=length(x), ncol=length(x))
  d = ncol(y)
  n   = sum(freq)
  m   = ncol(X)
  beta= matrix(param[1:(m*d)],m,d)
  D  = exp(-2*param[(m*d+1):(m*d+d)])
  A  = diag(d)
  i0 = m*d+d*(d+1)/2
  if(d>1) A[!lower.tri(A,diag=TRUE)] = param[(m*d+d+1):i0]
  eta   = param[(i0+1):(i0+d)]
  if(is.na(fixed.df)) df = exp(param[i0+d+1])
  else df = fixed.df
  Oinv  = t(A) %*% diag(D,d,d) %*% A
  u     = y-X %*% beta
  Q     = as.vector(apply((u %*% Oinv)*u,1,sum))
  L     = as.vector(u %*% eta)
  sf    = if(df<10000) sqrt((df+d)/(Q+df)) else sqrt((1+d/df)/(1+Q/df))
  t.    = L*sf
  dlogft= (-0.5)*(1+d/df)/(1+Q/df)
  dt.dL = sf
  dt.dQ = (-0.5)*L*sf/(Q+df)
  logT. = log.pt(t., df+d)
  dlogT.= exp(dt(t., df+d, log=TRUE)- logT.)
  u.freq= u*freq
  Dbeta = (-2* t(X) %*% (u.freq*dlogft) %*% Oinv
            - outer(as.vector(t(X) %*% (dlogT. * dt.dL* freq)), eta)
            - 2* t(X) %*% (dlogT.* dt.dQ * u.freq) %*% Oinv )
  Deta  = apply(dlogT.*sf*u.freq, 2, sum)
  if(d>1){
    M  = 2*( diag(D,d,d) %*% A %*% t(u * dlogft
                                      + u * dlogT. * dt.dQ) %*% u.freq)
    DA = M[!lower.tri(M,diag=TRUE)]
  }
  else DA= NULL
  M     = ( A %*% t(u*dlogft + u*dlogT.*dt.dQ) %*% u.freq %*% t(A))
  if(d>1) DD = diag(M) + 0.5*n/D
  else DD = as.vector(M + 0.5*n/D)
  grad = (-2)*c(Dbeta,DD*(-2*D),DA,Deta)
  if(is.na(fixed.df)) {
    df0= if(df<Inf) df else 1e8
    if(df0<10000){
      diff.digamma=  digamma((df0+d)/2) - digamma(df0/2)
      log1Q= log(1+Q/df0)
    }
    else
    {
      diff.digamma= log1p(d/df0)
      log1Q = log1p(Q/df0)
    }
    dlogft.ddf = 0.5 * (diff.digamma - d/df0
                         + (1+d/df0)*Q/((1+Q/df0)*df0) - log1Q)
    eps   = 1.0e-4
    df1 = df0 + eps
    sf1 = if(df0<10000) sqrt((df1+d)/(Q+df1)) else sqrt((1+d/df1)/(1+Q/df1))
    logT.eps = log.pt(L*sf1, df1+d)
    dlogT.ddf = (logT.eps-logT.)/eps
    Ddf   = sum((dlogft.ddf + dlogT.ddf)*freq)
    grad = c(grad, -2*Ddf*df0)
  }
  if(trace) cat("mst.dev.grad: norm is ",sqrt(sum(grad^2)),"\n")
  return(grad)
}
#-------------

st.2logL.profile=function(X=matrix(rep(1,n)), y, freq, trace=FALSE,
                           fixed.comp = c(ncol(X)+2, ncol(X)+3),
                           fixed.values = cbind(c(-4,4), log(c(1,25))),
                           npts=51/length(fixed.comp), plot.it=TRUE, ...)
{# plot2D profile deviance (=2(max.logL-logL)) using either parameters
  # if(plot.it & !exists(.Device)) stop("Device not active")
  #
  if(missing(freq)) freq = rep(1,length(y))
  n = sum(freq)
  m = ncol(X)
  npts = as.integer(npts)
  if(length(fixed.comp) == 1){
    param1 = seq(fixed.values[1], fixed.values[2], length=npts)
    logL = param2 = rep(NA,npts)}
  else{
    param1 = seq(fixed.values[1,1], fixed.values[2,1], length=npts)
    param2 = seq(fixed.values[1,2], fixed.values[2,2], length=npts)
    logL   = matrix(NA,npts,npts)}
  ls = lm.fit(X,y)
  omega = sqrt(var(resid(ls)))
  param = c(coef(ls), log(omega), 0, log(20))[-fixed.comp]
  max.logL = (-Inf)
  if(trace) cat(c("Running up to",npts,":"))
  for(i in 1:npts){
    if(trace) cat(" ",i)
    if(length(fixed.comp) == 1) {
      opt  = optim(param, fn=st.dev.fixed, method="Nelder-Mead",
                    X=X, y=y, freq=freq, trace=trace,
                    fixed.comp=fixed.comp, fixed.values=param1[i])
      logL[i] = opt$value/(-2)
      param = opt$par
      if(logL[i] > max.logL) {
        max.logL= logL[i]
        param = numeric(m+3)
        param[fixed.comp]  = param1[i]
        param[-fixed.comp] = opt$par
        dp= c(param[1:m], exp(param[m+1]), param[m+2], exp(param[m+3]))
        best = list(fixed.comp1=param1[i], fixed.comp2=NA,
                     dp=dp, logL=max.logL, opt=opt)
        param = param[-fixed.comp]
      }}
    else{
      for(j in 1:npts){
        opt  = optim(param, fn=st.dev.fixed, method="Nelder-Mead",
                      X=X, y=y, freq=freq, trace=trace,
                      fixed.comp=fixed.comp,
                      fixed.values=c(param1[i], param2[j] ))
        logL[i,j] = opt$value/(-2)
        if(j==1) param0 = opt$par
        if(j<npts)
          param = opt$par
        else
          param = param0
        if(logL[i,j] > max.logL) {
          max.logL= logL[i,j]
          param = numeric(m+3)
          param[fixed.comp]  = c(param1[i], param2[j])
          param[-fixed.comp] = opt$par
          dp= c(param[1:m], exp(param[m+1]), param[m+2], exp(param[m+3]))
          best = list(fixed.comp1=param1[i], fixed.comp2=param2[j],
                       dp=dp, logL=max.logL, opt=opt)
          param = param[-fixed.comp]
        }
      }}
  }
  if(trace) cat("\n")
  dev = 2 * (max(logL) - logL)
  if(plot.it){
    if(length(fixed.comp) == 1){
      plot(param1, dev, type="l", ...)
      points(x=best$fixed.comp1, y=0, pch=1)
    }
    else{
      contour(param1, param2, dev, labcex=0.5,
              levels=c(0.57, 1.37, 2.77, 4.6, 5.99, 9.2),
              labels=c(0.25, 0.5,  0.75, 0.90,0.95, 0.99),
              ...)
      points(x=best$fixed.comp1, y=best$fixed.comp2, pch=1,cex=0.5)
    }
  }
  title(main=paste("Dataset:", deparse(substitute(y)),
                   "\nProfile deviance", sep= " "))
  invisible(list(call=match.call(), param1=param1, param2=param2,
                 deviance=dev, max.logL=max.logL, best=best))
}


st.dev.fixed = function(free.param, X, y, freq, trace=FALSE,
                         fixed.comp=NA,  fixed.values=NA)
{# param components: beta, log(omega), alpha, log(df)
  n = sum(freq)
  m = ncol(X)
  param = numeric(length(free.param)+length(fixed.comp))
  param[fixed.comp]  = fixed.values
  param[-fixed.comp] = free.param
  beta  = param[1:m]
  omega = exp(param[m+1])
  eta = param[m+2]/omega
  df  = exp(param[m+3])
  u =  y - X %*% beta
  Q = freq*(u/omega)^2
  L = u*eta
  logDet = 2*log(omega)
  if(df < 10000)  {
    const= lgamma((df + 1)/2)- lgamma(df/2)-0.5*logb(df)
    log1Q = logb(1+Q/df)
  }
  else {
    const = (-0.5*logb(2)+ log1p((1/2)*(-1/2)/df))
    log1Q = log1p(Q/df)
  }
  dev = (n*(logDet - 2*const+ logb(pi)) + (df+1) * sum(freq * log1Q)
          -2*sum(log(2)+log.pt(L * sqrt((df+1)/(Q+df)),df+1)))
  if(trace) cat("st.dev.fixed (param, dev): ", param, dev,"\n")
  dev
}
log.pt = function(x, df){
  # fix for log(pt(...)) when it gives -Inf
  # see Abramowitz & Stegun formulae 26.7.8 & 26.2.13)
  # However, new releases of R (>=2.3) seem to have fixed the problem
  if(df == Inf) return(pnorm(x, log.p=TRUE))
  p = pt(x, df=df, log.p=TRUE)
  ninf = (p == -Inf)
  x0 = (1-1/(4*df))*(-x[ninf])/sqrt(1+x[ninf]^2/(2*df))
  p[ninf] = dnorm(x0,log=TRUE)-log(x0)+log1p(-1/(x0^2+2))
  p
}

num.deriv2 = function(x, FUN, ...)
{# derivate seconde numeriche, se FUN calcola il gradiente
  FUN = get(FUN, inherits = TRUE)
  values = FUN(x, ...)
  p = length(values)
  H = matrix(0, p, p)
  delta = cbind((abs(x) + 1e-10) * 1e-5, rep(1e-06, p))
  delta = apply(delta, 1, max)
  for(i in 1:p) {
    x1 = x
    x1[i] = x1[i]+delta[i]
    H[, i] = (FUN(x1, ...) - values)/delta[i]
  }
  (H+t(H))/2
}


msn.dev=function(param, X, y, freq, trace=FALSE)
{
  d = ncol(y)
  if(missing(freq)) freq=rep(1,nrow(y))
  n = sum(freq)
  m = ncol(X)
  beta=matrix(param[1:(m*d)],m,d)
  al.om=param[(m*d+1):(m*d+d)]
  y0 = y-X %*% beta
  Omega = (t(y0) %*% (y0*freq))/n
  D = diag(qr(2*pi*Omega)[[1]])
  logDet = sum(log(abs(D)))
  dev = n*logDet-2*sum(zeta(0,y0 %*% al.om)*freq)+n*d
  if(trace) {
    cat("\nmsn.dev:",dev,"\n","parameters:");
    print(rbind(beta,al.om))
  }
  dev
}
msn.marginal = function(xi=rep(0,length(alpha)), Omega, alpha,
                         comp=1:d, dp=NULL)
  
{# calcola parametri della marginale associata a comp di un SN_d
  # cfr SJS 2003, p.131-2
  if(!is.null(dp)){
    if(!is.null(dp$xi)) xi = dp$xi
    else
      if(!is.null(dp$beta)) xi = as.vector(dp$beta)
      Omega = dp$Omega
      alpha = dp$alpha
  }
  alpha = as.vector(alpha)
  d = length(alpha)
  xi = as.vector(xi)
  comp = as.integer(comp)
  if(length(xi) != d) stop("parameter size not compatible")
  if(all(dim(Omega) != c(d,d))) stop("parameter size not compatible")
  if(length(comp)<d){
    if(any(comp>d | comp<1)) stop("comp makes no sense")
    O   = cov2cor(Omega)
    O11 = O[comp,comp, drop=FALSE]
    O12 = O[comp,-comp, drop=FALSE]
    O21 = O[-comp,comp, drop=FALSE]
    O22 = O[-comp,-comp, drop=FALSE]
    alpha1 = matrix(alpha[comp], ncol=1)
    alpha2 = matrix(alpha[-comp], ncol=1)
    O11_inv = solvePD(O11)
    O22.1 = O22 - O21 %*% O11_inv %*% O12
    a.sum = as.vector(t(alpha2) %*% O22.1 %*% alpha2)
    a.new = as.vector(alpha1 + O11_inv %*% O12 %*% alpha2)/sqrt(1+a.sum)
    result= list(xi=xi[comp], Omega=Omega[comp,comp], alpha=a.new)
  }
  else {
    if(any(sort(comp)!=(1:d))) stop("comp makes no sense")
    result = list(xi=xi[comp], Omega=Omega[comp,comp], alpha=alpha[comp])
  }
  if(!is.null(dp$tau)) result$tau = dp$tau
  result
}




################################################################
## Below functions are only used in simulation studies


################################################################
#### data generation
## This function simulated bivariate data from a given distribution
## Inputs: real: real values of parameters
##         n: sample size
##         N: number of simulations
##         dist: name of distribution
################################################################
data_generate = function(real, n, N, dist, seed = 2018){
  
  set.seed(seed)
  X = matrix(0,nrow = n*2, ncol = N) 
  ## each column stands for a vector of simulated samples with the 1:n to be x and n:n+1 to be y
  
  if(dist == "log" || dist == "hr"){
    theta = real
    
    for(i in 1:N){
      X[,i] = as.vector(rbvevd(n, dep = theta, model = dist, mar1 = c(1,1,1)))}
  }
  
  if(dist == "alog"){
    theta = real[1]
    psi = real[2:3]
    
    for(i in 1:N){
      X[,i] = as.vector(rbvevd(n = n, dep = theta, 
                                asy = psi, model = dist, mar1 = c(1,1,1)))}
  }
  
  if(dist == "t"){
    nu = real[1]
    rho = real[2]
    
    for(i in 1:N){
      scale = matrix(c(1,rho,rho,1),2,2)
      X[,i] = as.vector(rmvt(n,sigma = scale, nu))}
  }

  return(X)
}


################################################################
#### main code to compute real values
## This function calculate values of eta_0, eta_star, VaR_Y, 
##                                   CoVaR_{Y|X}, CoVaR_{Y|X}_star
## Inputs: par: true values of parameter
##         family: parameteric model
##         p: risk level vector
##         VaR_X: conditional event in CoVaR
################################################################
real_compute = function(par, family, p){
  
  eta_star = eta_estimate(par_hat = par, p = p, family = family)
  
  if(family %in% c("log", "alog", "hr")){
    
    VaR_Y = qfrechet(1-p[2], loc=0, scale=1, shape=1, lower.tail = TRUE)
    
    CoVaR_star = qfrechet(1-p[2]*eta_star, loc=0, scale=1, shape=1, lower.tail = TRUE)
    
    Var_X = qfrechet(1-p[1], loc=0, scale=1, shape=1, lower.tail = TRUE)
    
    CoVaR = uniroot(covar_root, par = par, family = family, p = p, VaR_X = Var_X, 
                     interval = c(0, 1000000))$root
    
    eta = pfrechet(CoVaR, loc=0, scale=1, shape=1, lower.tail = FALSE)/p[2]
    }
  
  if(family %in% c("t")){
    
    VaR_Y = qt(1-p[2], df = par[1], lower.tail = TRUE)
    
    CoVaR_star = qt(1-p[2]*eta_star, df = par[1], lower.tail = TRUE)
    
    Var_X = qt(1-p[1], df = par[1], lower.tail = TRUE)
    
    CoVaR = uniroot(covar_root, par = par, family = family, p = p, VaR_X = Var_X, 
                     interval = c(0, 1000000))$root
    
    eta = pt(CoVaR, df = par[1], lower.tail = FALSE)/p[2]}
  
  return(c(eta = eta, eta_star = eta_star, VaR_Y = VaR_Y, CoVaR = CoVaR, CoVaR_star = CoVaR_star))
  
}

covar_root = function(y, par, family, p, VaR_X){
  
  
  if(family %in% c("log", "hr")){
    
    re = pbvevd(q = c(VaR_X, y), dep = par, model = family,
                 mar1 = c(1,1,1), lower.tail = FALSE) - p[1]*p[2]
    return(re)
  }
  
  if(family %in% c("alog")){
    
    re = pbvevd(q = c(VaR_X, y), dep =par[1], asy = par[2:3], model = family,
                 mar1 = c(1,1,1), lower.tail = FALSE) - p[1]*p[2]
    return(re)
  }
  
  if(family == "t"){
    
    scale = matrix(c(1, par[2], par[2],1),2,2)
    
    re = pmvt(lower = c(VaR_X, y), df = par[1], sigma = scale) - p[1]*p[2]
    
    return(re)
  }
}



################################################################
#### main code to test MRV

################################################################
## This function to test radial component
## Inputs: data: multivariate data frame
##         k: sample fraction (can in percentage)
##         eta: check Einmahl et al. [2021]
##         scale: scale the data before testing or not
##         alpha: significance level
################################################################
test_R = function(data, k = 0.1, eta = 0.5, scale = T, alpha){
  if(scale) data = apply(data, MARGIN = 2, FUN = scale)
  if(is.null(ncol(data))) R = data
  if(!is.null(ncol(data))) R = sqrt(rowSums(data^2))

  n = length(R)
  if(k < 1) k = round(k*n)
  R_k = sort(R, decreasing = F)[n-k]
  ind = which(R > R_k)
  gamma_all = sum(log(R[ind]) - log(R_k))/k
  int_t = function(t){
    R_kt = sort(R, decreasing = F)[n-round(k*t)]
    ((log(R_kt) - log(R_k))/gamma_all + log(t))^2*t^eta
  }
  Qn = k*integrate(int_t, lower = 0, upper = 1, subdivisions = 1000)$value
  if(alpha == 0.975) reject = ifelse(Qn > 1.792, 1, 0) ## 0.975
  if(alpha == 0.95) reject = ifelse(Qn > 1.515, 1, 0) ## 0.95
  if(alpha == 0.99) reject = ifelse(Qn > 2.145, 1, 0) ## 0.99
  return(c(Qn, reject))
} 

################################################################
## The functions to test angular component
## Inputs: data: multivariate data frame
##         df: dataframe to subset
##         colnum: which column to look at
##         nsplit: number of subsets
##         final_subsets: list used to store final subsets
##         k: sample fraction (can in percentage)
##         m: the number of cutpoint at each dim
##         scale: scale the data before testing or not
################################################################
subset_bycol = function(df, colnum, nsplit, final_subsets){
  # split the dataset into 3 pieces by specific column values
  #
  # Args:
  #   df: dataframe
  #     dataframe to subset
  #   colnum: integer
  #     which column to look at
  #   nsplit: integer
  #     number of subsets
  #   final_subsets: list
  #     list used to store final subsets
  # Return: list
  #   a list of n subsets
  break_points = quantile(df[,colnum], seq(0,1,1/nsplit))
  res = list()
  for (i in 1:nsplit) {
    if (i == 1) {
      subset_row = which(df[,colnum] == break_points[i+1])
      subset = df[subset_row,]
      res[[i]] = subset
      
    } else {
      subset_row = which((df[,colnum] > break_points[i]) & (df[,colnum] = break_points[i+1]))
      subset = df[subset_row,]
      res[[i]] = subset
    }
  }
  if (colnum == (ncol(df) -1)) { # splitted at the last column
    for (sub_df in res) {
      final_subsets[[length(final_subsets)+1L]] = sub_df
    }
    return(final_subsets)
  }
  for (i in 1:nsplit) { # recursion
    final_subsets = Recall(res[[i]], colnum+1, nsplit, final_subsets)
  }
  return(final_subsets)
}

subset_main = function(df, nsplit) {
  # split a dataframe into nsplit^nvar subsets
  #
  # Arg:
  #   df: dataframe to split
  #   nsplit: number of splits for each column
  # 
  # Return:
  #   a list of nsplit^nvar subsets
  res = list()
  final_res = subset_bycol(df, 1, nsplit, res)
  return(final_res)
}

test_A = function(data, k = 0.1, m, scale = F){
  # m is the number of cutpoint, the k-largest points are divided into m disjoint parts with about equal number of obs.
  dim = ncol(data)
  if(scale) data = apply(data, MARGIN = 2, FUN = scale)
  R = sqrt(rowSums(data^2))
  Theta = atan2(data[,2], data[,1]) 
  #indd = NULL                 ## split 1
  #indd = which(Theta < -pi/2) ## split 2
  #indd = which(Theta < 0)     ## split 3
  #ind = which(Theta < pi/2)   ## split 3
  #Theta[indd] = Theta[indd] + 2*pi
  n = length(R)
  if(k < 1) k = round(k*n)
  R_k = sort(R, decreasing = F)[n-k]
  ind = which(R > R_k)
  gamma_all = sum(log(R[ind]) - log(R_k))/k
  
  sub = data.frame(cbind(Theta[ind], R[ind]))
  re = subset_main(df = sub, nsplit = m)
  gamma = rep(0, length(re))
  for(i in 1:length(re)){
    gamma[i] = sum(log(re[[i]][,ncol(sub)]) - log(R_k))/nrow(re[[i]])
  }
  Tn = (k/m^(dim-1))*sum((gamma/gamma_all - 1)^2)
  pvalue = pchisq(Tn, df = (m^(dim-1)-1), lower.tail = F)
  reject = ifelse(pvalue < 0.025, 1,0)
  return(list(subset = re, test = data.frame(Tn = Tn, pvalue = pvalue, reject = reject)))
}