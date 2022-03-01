##############################################################
############## Dynamic forecast ########################
path = ".../"
source(paste(path, "Functions_CoVaR.R", sep = ""))
##############################################################

###### Loading data ############################
Dat = fread(paste(path,"Snp500_dataset.csv",sep=""))
Dat$date <- as.Date(Dat$date)
Ldata <- data.table(dcast(Dat,formula=date ~ symbol, value.var="price", drop = FALSE))
Ldata$date <- as.Date(Ldata$date)
tsdata = as.timeSeries(Ldata)
Losses = Return.calculate(tsdata, method = "log")[-1,] * (-100)
Date = Dat$date[which(Dat$symbol == "^GSPC")][-1]
rm(Dat, tsdata)

Losses = data.frame(Losses)
colnames(Losses) <- c(colnames(Losses)[-15], "GSPC")

####### Set parametes ##########################
level <- c(0.02,0.05)
length.out <- nrow(Losses) - 3000
begin <- round(seq(1, length.out, 50))

m_group <-  round(c(0.09, 0.14, 0.09, 0.08, 0.03)*3000)
dist_group <- c("log", "hr", "bilog", "alog", "t")
k_sys <- c(110, 110)
k_ins <- c(150, 120)

#########CoVaR estimation#######################

ins <- 1 ## start with AFL
SData <- cbind(Losses[,ins], Losses$Nasdaq)
CoVaR_est <- matrix(0, nrow = 50*length(begin), ncol = (length(dist_group) + 2))
VaR_ins_est <- rep(0, 50*length(begin))
for(i in 1:length(begin)){
  start <- begin[i]
  end <- begin[i] + 2999 + 50
  
  if(end > nrow(SData)) end <- nrow(SData)
  Sfilter_ins <- filtering(SData[start:end,1], forecast = TRUE, n_out = 50)
  Sfilter_sys <- filtering(SData[start:end,2], forecast = TRUE, n_out = 50)
  
  Sinnos <- cbind(Sfilter_ins$residuals, Sfilter_sys$residuals)
  
  CoVaR_Z <- covar_est(Sinnos, group = dist_group, m = m_group, k = k_sys, p = level)## proposed method
  
  ## other two methods
  VaR_ins_Z <- VaR_estimate(dat = Sinnos[,1], k = k_ins, p = level[1])
  CoVaR_other_Z <- CoVaR_other(Data = Sinnos, VaR = VaR_ins_Z, level = level)
  
  VaR_ins_est[(50*i-49):(50*i)] <- VaR_ins_Z*Sfilter_ins$sigma.forecast + Sfilter_ins$mean.forecast
  
  CoVaR_est[(50*i-49):(50*i),] <- matrix(rep(c(CoVaR_Z, CoVaR_other_Z), 50), nrow = 50, byrow = T)*Sfilter_sys$sigma.forecast 
  + Sfilter_sys$mean.forecast
  
  VaR_ins_est <- VaR_ins_est[1:length.out]
  CoVaR_est <- CoVaR_est[1:length.out,]
}
colnames(CoVaR_est) <- c("log", "hr", "bilog", "alog", "t", "FP", "EVT_NZ")

#### Average score #############################
test <- SData[(length.fit + 1):nrow(SData),]
average_scores <- average_score(CoVaR_est, VaR_ins_est, test = test, level = level[2])


### Traffic light matrix ######################
## repeat CoVaR estimation for all firms to get all estimates
## following similar steps as in insample analysis