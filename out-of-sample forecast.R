##############################################################
############## Dynamic forecast ########################
path = ".../"
source(paste(path, "Functions_CoVaR.R", sep = ""))
##############################################################
library(data.table)
library(timeSeries)
library(PerformanceAnalytics)

###### Loading data ############################
Dat = fread(paste(path,"Snp500_dataset.csv",sep=""))
Dat$date = as.Date(Dat$date)
Ldata = data.table(dcast(Dat,formula=date ~ symbol, value.var="price", drop = FALSE))
Ldata$date = as.Date(Ldata$date)
tsdata = as.timeSeries(Ldata)
Losses = Return.calculate(tsdata, method = "log")[-1,] * (-100)
Losses = data.frame(Losses)
rownames(Losses) = Ldata$date[-1]
Losses = data.frame(Losses[,-1], GSPC = Losses[,1])
rm(Dat, tsdata, Ldata)


####### Set parameters ##########################
level = c(0.02,0.05)
length.out = nrow(Losses) - 3000
begin = round(seq(1, length.out, 50))

m_group =  round(c(0.09, 0.14, 0.09, 0.08, 0.03)*3000)
dist_group = c("log", "hr", "bilog", "alog", "t")
k_sys = ## a sequence of best values selected by minimizing in-sample ave score

#########CoVaR estimation#######################

ins = 1 ## start with AFL
k_ins = ## a sequence of best values selected by minimizing in-sample ave score
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
  
  CoVaR_Z = covar_est(Sinnos, group = dist_group, m = m_group, k = c(150,k_sys[j]), 
                       p = level)## proposed method
  
  ## other two methods
  VaR_ins_Z = VaR_estimate(dat = Sinnos[,1], k = c(150,k_ins[j]), p = level[1])
  CoVaR_NZ_Z = CoVaR_NZ(Data = Sinnos, VaR = VaR_ins_Z, level = level)
  CoVaR_FP_Z = CoVaR_NZ(Data = Sinnos, fit_par = Sfilter_ins$coef, level = level)
  
  VaR_ins_est[(50*i-49):(50*i)] = VaR_ins_Z*Sfilter_ins$sigma.forecast + Sfilter_ins$mean.forecast
  
  CoVaR_est[(50*i-49):(50*i),] = matrix(rep(c(CoVaR_Z, CoVaR_FP_Z, CoVaR_NZ_Z), 50), nrow = 50, byrow = T)*Sfilter_sys$sigma.forecast 
  + Sfilter_sys$mean.forecast
  
  VaR_ins_est = VaR_ins_est[1:length.out]
  CoVaR_est = CoVaR_est[1:length.out,]
}
colnames(CoVaR_est) = c("log", "hr", "bilog", "alog", "t", "FP", "EVT_NZ")

#### Average score #############################
test = SData[(length.fit + 1):nrow(SData),]
average_scores = average_score(CoVaR_est, VaR_ins_est, test = test, level = level[2])
colnames(ave_score) = c("log", "hr", "bilog", "alog", "t", "FP", "EVT_NZ")
rownames(ave_score) = colnames(Losses)[-16]

### Traffic light matrix ######################
###### repeat above steps for each institution
###### we can store the CoVaR estimates into a list call CoVaR_output; each component stands for a institution
###### and store the VaR estimates into a matrix call VaR_output; each column stands for a institution

####### Test with all institutions

pick_out <- function(CoVaR_list, i){
  ## pick estimates of CoVaR of i-th method out for all institutions
  re <- CoVaR_list[[1]][,i]
  for(s in 2:length(CoVaR_list)){
    re = cbind(re, CoVaR_list[[s]][,i])
  }
  return(re)
}

M_01 <- diag(rep(NA, 7))
colnames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
rownames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
M_05 <- M_10 <- Tn <- M_01
for(i in 1:7){
  for(j in 1:7){
    if(i == j) M_01[i,j] = NA
    if(i!= j){
      re = comparative.test(VaR_fr = VaR_output, CoVaR_fr = pick_out(CoVaR_output, i), 
                            CoVaR_fr_r = pick_out(CoVaR_output, j), 
                            test_ins = Losses[3001:5534,-ncol(Losses)], 
                            test_sys = Losses[3001:5534,"GSPC"], p = level[2])
      M_01[i,j] = re[1]
      M_05[i,j] = re[2]
      M_10[i,j] = re[3]
      
      Tn[i,j] = re[4]
    }
  }
}

melted_cormat <- melt(M_10, na.rm = TRUE)

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

Ind_comparative.test <- function(VaR_fr, CoVaR_fr, CoVaR_fr_r, 
                                 test_ins, test_sys){
  ind_pick <- which(test_ins > VaR_fr)
  test_sub <- test_sys[ind_pick]
  CoVaR_fr_sub <- CoVaR_fr[ind_pick]
  CoVaR_fr_r_sub <- CoVaR_fr_r[ind_pick]
  tmp = tmpS = list()
  tmp[[1]] <- score_fun(x = test_sub, r = CoVaR_fr_sub, level = 0.05)
  tmpS[[1]] <- score_fun(x = test_sub, r = CoVaR_fr_r_sub, level = 0.05)
  re <- TLMfn.biv(tmp, tmpS, cor = TRUE)
  return(re)
}
Compute_M <- function(ins, plot = 10){
  M_01 <- diag(rep(NA, 7))
  colnames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
  rownames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
  M_05 <- M_10 <- Tn <- M_01
  for(i in 1:7){
    for(j in 1:7){
      if(i == j) M_01[i,j] = NA
      if(i!= j){
        re = Ind_comparative.test(VaR_fr = VaR_output[,ins], CoVaR_fr = CoVaR_output[[ins]][,i], 
                                  CoVaR_fr_r = CoVaR_output[[ins]][,j], 
                                  test_ins = Losses[3001:5534,ins], 
                                  test_sys = Losses[3001:5534,"^GSPC"])
        M_01[i,j] = re[1]
        M_05[i,j] = re[2]
        M_10[i,j] = re[3]
        
        Tn[i,j] = re[4]
      }
    }
  }
  if(plot == 5) melted_cormat <- melt(M_05, na.rm = TRUE)
  if(plot == 1) melted_cormat <- melt(M_01, na.rm = TRUE)
  if(plot == 10) melted_cormat <- melt(M_10, na.rm = TRUE)
  
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