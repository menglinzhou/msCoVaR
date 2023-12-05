##############################################################
############## Model diagnostic ########################

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

###### Get realized innovations ################
innovation = Losses
for(i in 1:ncol(Losses)){
  innovation[,i] = filtering(Losses[,i], forecast = FALSE)
}

###### acf plot ###############################
for(i in 1:length(select)){
  acf(innovations[,i], ylim = c(-0.05,0.1), main = select[i])
}
for(i in 1:length(select)){
  acf((innovations[,i])^2, ylim = c(-0.05,0.1),main = select[i])
}

###### Ljung-Box test  #######################
pvalue20 = pvalue30 = rep(NA, length(select))
for(i in 1:length(select)){
  re = Box.test(innovations[,i],lag=20,type = "Ljung-Box")
  pvalue20[i] = re$p.value
  re = Box.test(innovations[,i],lag=30,type = "Ljung-Box")
  pvalue30[i] = re$p.value
}

##### Estimate upper tail coefficients #####
lambda_est = rep(NA, ncol(innovation)-1)
for(i in 1:(ncol(innovation)-1)){
  lambda_est[i] = tailDep(u1 = innovations[,i], u2 = innovation$GSPC, 
                          alpha = 10:30, 
                          rank = FALSE, eps = 0.2,iprint = F)
}

##### Estimate tail index #################
# k is selected from sensitive analysis
k = c(250,250,300,250,280,250,250,250,270,300,250,250,250,250,250,200)
gamma_est = rep(NA, ncol(innovation))
for(i in 1:length(gamma_est)){
  gamma_est[i] = tail_estimate(dat = innovation[,i], k = k[i])
}


##### Test MRV ############################
k_seq = seq(50, 500, 10)
test_re = matrix(NA, nrow = length(k_seq), ncol = 2)
for(i in 1:length(k_seq)){
  test_re[i,] = test_R(data = innovations$GSPC, k = k_seq[i], scale = FALSE, 
                       alpha = 0.95)
}
test_re = cbind(k_seq, test_re)
colnames(test_re) <- c("k", "test stat", "reject")

