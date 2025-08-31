##############################################################
############## Dynamic forecast ########################
path = getwd()
source(paste(path, "/code/functions.R", sep = ""))
##############################################################
library(data.table)
library(timeSeries)
library(PerformanceAnalytics)

###### Loading data ############################
Dat = fread(paste(path,"/data/Snp500_dataset.csv",sep=""))
Dat$date = as.Date(Dat$date)
Ldata = data.table(dcast(Dat,formula=date ~ symbol, value.var="price", drop = FALSE))
Ldata$date = as.Date(Ldata$date)
tsdata = as.timeSeries(Ldata)
Losses = Return.calculate(tsdata, method = "log")[-1,] * (-100)
Losses = data.frame(Losses)
rownames(Losses) = Ldata$date[-1]
Losses = data.frame(Losses[,-ncol(Losses)], GSPC = Losses[,ncol(Losses)])
rm(Dat, tsdata, Ldata)

####### Set parameters ##########################
level = c(0.02,0.05)
length.out = nrow(Losses) - 3000
begin = round(seq(1, length.out, 50))

m_group =  round(c(0.09, 0.14, 0.08, 0.03)*3000)
dist_group = c("log", "hr", "alog", "t")
k1 = 150
k2_sys = ## a sequence of best values selected by minimizing in-sample ave score (for proposed method)
###################################### CoVaR estimation with given parameters

ins = 1 ## start with AFL
k2_ins = ## a sequence of best values selected by minimizing in-sample ave score (for NZ method)
SData = cbind(Losses[,ins], Losses$GSPC)
CoVaR_est = matrix(0, nrow = 50*length(begin), ncol = (length(dist_group) + 2))
VaR_ins_est = rep(0, 50*length(begin))
for(i in 1:length(begin)){
  start = begin[i]
  end = begin[i] + 2999 + 50
  
  if(end > nrow(SData)) end = nrow(SData)
  Sfilter_ins = filtering(SData[start:end,1], forecast = TRUE, n_out = 50)
  Sfilter_sys = filtering(SData[start:end,2], forecast = TRUE, n_out = 50)
  
  Sinnos = cbind(Sfilter_ins$residuals, Sfilter_sys$residuals)
  
  CoVaR_Z = CoVaR_est(Sinnos, group = dist_group, m = m_group, 
                      k = c(k1,k2_sys[i]), 
                       p = level) ## proposed method
  
  ## other two methods
  tail_ins_Z = tail_estimate(dat = Sinnos[,1], k=k1)
  VaR_ins_Z = Qvar_esti(dat = Sinnos[,1], gamma = tail_ins_Z, k = k2_ins[i], 
                        p = level[1])
  CoVaR_NZ_Z = CoVaR_NZ(Data = Sinnos, VaR = VaR_ins_Z, level = level)
  CoVaR_FP_Z = CoVaR_FP(Data = Sinnos, fit_par = Sfilter_ins$coef, level = level)
  
  VaR_ins_est[(50*i-49):(50*i)] = VaR_ins_Z*Sfilter_ins$sigma.forecast + Sfilter_ins$mean.forecast
  
  CoVaR_est[(50*i-49):(50*i),] = matrix(rep(c(CoVaR_Z, CoVaR_FP_Z, CoVaR_NZ_Z), 50), nrow = 50, byrow = T)*Sfilter_sys$sigma.forecast + Sfilter_sys$mean.forecast
  
  VaR_ins_est = VaR_ins_est[1:length.out]
  CoVaR_est = CoVaR_est[1:length.out,]
}
colnames(CoVaR_est) = c("log", "hr", "alog", "t", "FP", "EVT_NZ")

#### Unconditional test and Average score #####
test = SData[(length.fit + 1):nrow(SData),]
test_result = unconditional_test(CoVaR_est, VaR_ins_est, test = test, level = level)
average_scores = average_score(CoVaR_est, VaR_ins_est, test = test, level = level[2])
colnames(ave_score) = c("log", "hr", "alog", "t", "FP", "EVT_NZ")
rownames(ave_score) = colnames(Losses)[-ncol(Losses)]

### Traffic light matrix ######################
###### repeat above steps for each institution
###### we can store the CoVaR estimates into a list call CoVaR_output; each component stands for a institution
###### and store the VaR estimates into a matrix call VaR_output; each column stands for a institution

####### Test with all institutions

pick_out = function(CoVaR_list, i){
  ## pick estimates of CoVaR of i-th method out for all institutions
  re = CoVaR_list[[1]][,i]
  for(s in 2:length(CoVaR_list)){
    re = cbind(re, CoVaR_list[[s]][,i])
  }
  return(re)
}

M_01 = diag(rep(NA, 7))
colnames(M_01) = c("Log", "HR", "Alog","t", "FP", "EVT-NZ")
rownames(M_01) = c("Log", "HR", "Alog","t", "FP", "EVT-NZ")
M_05 = M_10 = Tn = M_01
for(i in 1:7){
  for(j in 1:7){
    if(i == j) M_01[i,j] = NA
    if(i!= j){
      re = comparative.test(VaR_fr = VaR_output, CoVaR_fr = pick_out(CoVaR_output, i), 
                            CoVaR_fr_r = pick_out(CoVaR_output, j), 
                            test_ins = Losses[3001:nrow(Losses),-ncol(Losses)], 
                            test_sys = Losses[3001:nrow(Losses),"GSPC"], p = level[2])
      M_01[i,j] = re[1]
      M_05[i,j] = re[2]
      M_10[i,j] = re[3]
      
      Tn[i,j] = re[4]
    }
  }
}

melted_cormat = melt(M_10, na.rm = TRUE)

ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "yellow", breaks = NULL,
                       midpoint = 1, limit = c(0,2), space = "Lab") +
  labs(title = "p-value = 0.1", x = "Reference method", y = "Proposed method") + 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 0, size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))


######## Test with a single institution

Ind_comparative.test = function(VaR_fr, CoVaR_fr, CoVaR_fr_r, 
                                 test_ins, test_sys){
  ind_pick = which(test_ins > VaR_fr)
  test_sub = test_sys[ind_pick]
  CoVaR_fr_sub = CoVaR_fr[ind_pick]
  CoVaR_fr_r_sub = CoVaR_fr_r[ind_pick]
  tmp = tmpS = list()
  tmp[[1]] = score_fun(x = test_sub, r = CoVaR_fr_sub, level = 0.05)
  tmpS[[1]] = score_fun(x = test_sub, r = CoVaR_fr_r_sub, level = 0.05)
  re = TLMfn.biv(tmp, tmpS, cor = TRUE)
  return(re)
}
Compute_M = function(ins, plot = 10){
  M_01 = diag(rep(NA, 7))
  colnames(M_01) = c("Log", "HR", "Alog","t", "FP", "EVT-NZ")
  rownames(M_01) = c("Log", "HR", "Alog","t", "FP", "EVT-NZ")
  M_05 = M_10 = Tn = M_01
  for(i in 1:7){
    for(j in 1:7){
      if(i == j) M_01[i,j] = NA
      if(i!= j){
        re = Ind_comparative.test(VaR_fr = VaR_output[,ins], CoVaR_fr = CoVaR_output[[ins]][,i], 
                                  CoVaR_fr_r = CoVaR_output[[ins]][,j], 
                                  test_ins = Losses[3001:5534,ins], 
                                  test_sys = Losses[3001:5534,"GSPC"])
        M_01[i,j] = re[1]
        M_05[i,j] = re[2]
        M_10[i,j] = re[3]
        
        Tn[i,j] = re[4]
      }
    }
  }
  if(plot == 5) melted_cormat = melt(M_05, na.rm = TRUE)
  if(plot == 1) melted_cormat = melt(M_01, na.rm = TRUE)
  if(plot == 10) melted_cormat = melt(M_10, na.rm = TRUE)
  
  p = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "red", high = "green", mid = "yellow", breaks = NULL,
                         midpoint = 1, limit = c(0,2), space = "Lab") +
    labs(title = ins, x = "Reference method", y = "Proposed method") + 
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 0, size = 13), 
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15))
  return(p)
}

Compute_M("UNM")


####### Below is the procedure to select "best" k2 for the proposed method
####### Same procedure is applied to NZ method to select "best" k2

####### Apply GARCH filtering ###################
filter_all = function(dat, model = "sGARCH"){
  ### this function applies GARCH filtering to all 51 samples with size 3000
  innos = matrix(0, nrow = 3000, ncol = length(begin))
  coef_est = matrix(0, ncol = 7, nrow = length(begin))
  mean_f = sigma_f = NA
  for(i in 1:length(begin)){
    start = begin[i]
    end = begin[i] + 2999 + 50
    if(end > length(dat)) end = length(dat)
    filter_re = tryCatch(filtering(dat[start:end], forecast = TRUE, n_out = (end-start-2999), 
                                    model = model),
                          error = function(e){
                            filtering(dat[start:end], forecast = TRUE, 
                                      n_out = (end-start-2999), 
                                      model = "sGARCH")})
    innos[,i] = filter_re$residuals
    coef_est[i,] = filter_re$coef
    mean_f = c(mean_f, filter_re$mean.forecast)
    sigma_f = c(sigma_f, filter_re$sigma.forecast)
    print(i)
  }
  mean_f = mean_f[-1]
  sigma_f = sigma_f[-1]
  return(list(innovations = innos, coef = coef_est,
              mean.forecasts = mean_f, sigma.forecasts = sigma_f))
}

cl = parallel::makeCluster(getOption('cl.cores', 6))
doParallel::registerDoParallel(cl)
filter_result = foreach(j=1:ncol(Losses), .export = ls(envir = globalenv()), 
                         .packages=c('rugarch'), .errorhandling = "pass") %dopar% 
  filter_all(Losses[,j])
stopCluster(cl)
names(filter_result) = colnames(Losses)

####### M-estimation ############################
m_group =  round(c(0.09, 0.14, 0.08, 0.03)*3000)
dist_group = c("log", "hr", "alog", "t")

M_est_allinnos = function(ins, sys, dist, m){
  ### This function do M-estimation for all 51 samples with size 3000 for an ins
  cl = parallel::makeCluster(getOption('cl.cores', 6))
  doParallel::registerDoParallel(cl)
  M_est = foreach(i = 1:ncol(ins), .export = ls(envir = globalenv()), 
                   .combine = rbind,
                   .packages=c('rugarch'), .errorhandling = "pass") %dopar% 
    M_estimate(X = cbind(ins[,i], sys[,i]), family = dist, m = m)
  stopCluster(cl)
  return(M_est)
}

M_result = list()
for(j in 1:num_ins){
  re = list()
  for(s in 1:length(dist_group)){
    re[[s]] = M_est_allinnos(ins = filter_result[[j]]$innovations, sys = filter_result[["GSPC"]]$innovations,
                             dist = dist_group[s], m = m_group[s])
    names(re) = dist_group
  }
  M_result[[j]] = re
  names(M_result) = names(filter_result)[1:j]
}
rm(re, j, s)


####### Estimate VaR_Y with different k2 ######
k1 = 150
tail_sys = rep(NA, length(begin))
for (i in 1:length(tail_sys)){
  tail_sys[i] = tail_estimate(filter_result[["GSPC"]]$innovations[,i], k = 150)
}


k2_seq = seq(200, 300, 10)
VaR_sys = matrix(NA, nrow = length(begin), ncol = length(k2_seq))
for (i in 1:nrow(VaR_sys)){
  for (j in 1:ncol(VaR_sys))
    VaR_sys[i, j] = Qvar_esti(dat = filter_result[["GSPC"]]$innovations[,i], 
                              gamma = tail_sys[i], 
                              k = k2_seq[j], p = level[2])
}

####### Estimate CoVaR with different VaR_Y ######
get_CoVaR = function(ins, sys, VaR_sys, tail_sys, M_re, 
                      dist_group, level){
  
  re = matrix(0, ncol = length(dist_group), nrow = ncol(ins))
  
  for(i in 1:ncol(ins)){
    
    eta_est = rep(0, length(dist_group))
    
    for(j in 1:length(dist_group)){
      
      eta_est[j] = eta_estimate(par_hat = as.numeric(M_re[[j]][i,]), 
                                 p = level, 
                                 family = dist_group[j])
      
    }
    
    re[i,] = VaR_sys[i]*(eta_est)^(-tail_sys[i])
    
  }
  colnames(re) = dist_group
  return(re)
}

get_CoVaR_seq = function(ins, sys, VaR, tail, M_est){
  Log = HR = Alog = t = matrix(NA, nrow = nrow(VaR), ncol = ncol(VaR))
  for(i in 1:ncol(VaR)){
    re = get_CoVaR(ins = ins, sys = sys, VaR_sys = VaR[i,], 
                   tail_sys = tail,
                   M_re = M_est, dist_group = dist_group, 
                   level = c(0.02,0.05))
    Log[,i] = re[,1]
    HR[,i] = re[,2]
    Alog[,i] = re[,3]
    t[,i] = re[,4]
  }
  return(list(Log = Log, HR = HR, Alog = Alog, t = t))
}

cl = parallel::makeCluster(getOption('cl.cores', 6))
doParallel::registerDoParallel(cl)
all_CoVaR = foreach(i=1:num_ins, .export = ls(envir = globalenv()), 
                     .packages=c('evir', 'Matrix', 'MASS', 'sn', 'cubature'), 
                     .errorhandling = "pass") %dopar% 
  get_CoVaR_seq(ins = filter_result[[i]]$innovations, 
                sys = filter_result[["GSPC"]]$innovations,
                VaR = VaR_sys, tail = tail_sys, 
                M_est = M_result[[i]])
stopCluster(cl)
names(all_CoVaR) = names(M_result)

####### Get score for CoVaR at different k2 ######

get_score = function(CoVaR_est, VaR_ins, mean.f.sys, se.f.sys,
                      mean.f.ins, se.f.ins, test){
  ave_score = rep(NA, ncol(CoVaR_est))
  if(is.null(ncol(VaR_ins))){
    VaR_org = rep(VaR_ins, each = 50)[1:length.out]*se.f.ins + mean.f.ins
  }
  for(i in 1:ncol(CoVaR_est)){
    CoVaR_org = rep(CoVaR_est[,i], each = 50)[1:length.out]*se.f.sys + mean.f.sys
    
    if(!is.null(ncol(VaR_ins))){
      VaR_org = rep(VaR_ins[,i], each = 50)[1:length.out]*se.f.ins + mean.f.ins
    }
    
    ave_score[i] = average_score(forecast_CoVaR = matrix(CoVaR_org, ncol = 1), 
                                 forecast_VaR = VaR_org, test = test, level = level[2])
  }
  return(ave_score)
}

ave_score = list()
for(i in 1:num_ins){
  ## for each institution
  ins_name = colnames(Losses)[i]
  score_re = list()
  for(j in 1:length(dist_group)){
    ## for each parametric model
    score_mat = matrix(NA, nrow = length(seq_k), ncol = length(seq_k))
    for(s in 1:length(seq_k)){
      score_mat[,s] = get_score(CoVaR_est = all_CoVaR[[i]][[j]], 
                                VaR_ins = all_tail_VaR[[i]]$VaR[,s],
                                mean.f.sys = filter_result[["GSPC"]]$mean.forecasts, 
                                se.f.sys = filter_result[["GSPC"]]$sigma.forecasts,
                                mean.f.ins = filter_result[[i]]$mean.forecasts,
                                se.f.ins = filter_result[[i]]$sigma.forecasts,
                                test = Losses[3001:nrow(Losses),c(ins_name,"GSPC")])
    }
    score_re[[j]] = score_mat
  }
  names(score_re) = dist_group
  ave_score[[i]] = score_re
}
names(ave_score) = select
rm(score_re, i, j, score_mat, s)

### we minimize above score to select "best" k2