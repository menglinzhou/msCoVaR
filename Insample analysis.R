##############################################################
################## Insample Analysis #########################
path = ".../"
source(paste(path, "Functions_CoVaR.R", sep = ""))

##############################################################

###### Loading data ################################
Dat = fread(paste(path,"top10_nas.csv",sep=""))
Ldata <- data.table(dcast(Dat,formula=date ~ symbol, value.var="price", drop = FALSE))
Ldata$date <- as.Date(Ldata$date)
tsdata = as.timeSeries(Ldata)
Losses = Return.calculate(tsdata, method = "log")[-1,] * (-100)
Date = Date = Dat$date[which(Dat$symbol == "NFLX")][-1]
rm(Dat, tsdata)

Losses = data.frame(Losses)
Losses = na.omit(Losses)
Losses = Losses[, c("NFLX", "BKNG", "ILMN", "ISRG", "AMZN", "INTU", "ADBE", "REGN", "CSCO", "ATVI", "Nasdaq")]

#### Data filtering ################################
innovations <- Losses
for(i in 1:ncol(Losses)){
  innovations[,i] <- filtering(Losses[,i], forecast = FALSE)
}
colnames(innovations) <- colnames(Losses)

#### Set parameters ###############################
level <- c(0.02,0.05)
m_group <-  round(c(0.09, 0.14, 0.09, 0.08, 0.03)*nrow(innovations))
dist_group <- c("log", "hr", "bilog", "alog", "t")
k_sys <- c(200, 250)
k_ins <- matrix(rep(c(200,250),10), byrow = T, nrow = 10)

#### Insample CoVaR estimation #####################
CoVaR_insample <- matrix(0, nrow = (ncol(Losses) - 1), ncol = length(dist_group))
CoVaR_insample_other <- matrix(0, nrow = (ncol(Losses) - 1), ncol = 2)
VaR_ins <- rep(0, (ncol(Losses) - 1))
for(i in 1:(ncol(Losses) - 1)){
  CoVaR_insample[i,] <- covar_est(Data = cbind(innovations[,i], innovations$Nasdaq),
                                  group = dist_group, m = m_group, k = k_sys, p = level)
  VaR_ins[i] <- VaR_estimate(dat = innovations[,i], k = k_ins[i,], p = level[1])
  CoVaR_insample_other[i,] <- CoVaR_other(Data = cbind(innovations[,i], innovations$Nasdaq),
                                          VaR = VaR_ins[i], level = level)
  print(i)
}

CoVaR_output <- cbind(CoVaR_insample, CoVaR_insample_other)
colnames(CoVaR_output) <- c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")

#### Unconditional test & Average score ############
test_output <-list()
ave_score <- matrix(0, nrow = nrow(CoVaR_output), ncol = (length(dist_group) + 2))
for(i in 1:nrow(CoVaR_output)){
  test <- cbind(innovations[,i], innovations$Nasdaq)
  forecast_CoVaR <- matrix(rep(CoVaR_output[i,], nrow(test)), nrow = nrow(test), byrow = T)
  forecast_VaR <- rep(VaR_ins[i], nrow(test))
  test_output[[i]] <- unconditional_test(forecast_CoVaR, forecast_VaR, test = test, level = level)
  ave_score[i,] <- average_score(forecast_CoVaR, forecast_VaR, test = test, level = level[2])
}

#### Traffic light matrix ##########################

## pooled comparative test
M_01 <- diag(rep(NA, 7))
colnames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
rownames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
M_05 <- M_10 <- Tn <- M_01
for(i in 1:7){
  for(j in 1:7){
    if(i == j) M_01[i,j] = NA
    if(i!= j){
      re = pool_comparative.test(VaR_fr = VaR_ins, CoVaR_fr = CoVaR_output[,i], 
                                 CoVaR_fr_r = CoVaR_output[,j], 
                                 test_ins = innovations[,-11], test_sys = innovations[,11])
      M_01[i,j] = re[1]
      M_05[i,j] = re[2]
      M_10[i,j] = re[3]
      
      Tn[i,j] = re[4]
    }
  }
}

## plot the matrix
library(reshape2)
library(ggplot2)

melted_cormat <- melt(M_05, na.rm = TRUE)

ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "yellow", breaks = NULL,
                       midpoint = 1, limit = c(0,2), space = "Lab") +
  labs(x = "Reference method", y = "Proposed method") + 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 0, size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
coord_fixed()


## individual comparative test
Compute_M <- function(ins, plot = 5){
  M_01 <- diag(rep(NA, 7))
  colnames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
  rownames(M_01) = c("Log", "HR", "Bilog", "Alog","t", "FP", "EVT-NZ")
  M_05 <- M_10 <- Tn <- M_01
  for(i in 1:7){
    for(j in 1:7){
      if(i == j) M_01[i,j] = NA
      if(i!= j){
        re = Ind_comparative.test(VaR_fr = VaR_est_hill[ins], CoVaR_fr = CoVaR_output[ins,i], 
                                  CoVaR_fr_r = CoVaR_output[ins,j], 
                                  test_ins = innovations[,ins], test_sys = innovations[,11])
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
    labs(x = "Reference method", y = "Proposed method") + 
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 0, size = 13), 
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15))
  coord_fixed()
  # return(list(M_01 = M_01, M_05 = M_05, M_10 = M_10))
  return(p)
}

par(mfrow = c(2,3))
for(i in 1:10){
  Compute_M(colnames(innovations)[i])
}