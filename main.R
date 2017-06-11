#setwd("D:/Dropbox/Ricardo Lindo/M2/Thesis")
setwd('~/Desktop')
library(glmnet)
library(parallel)
# Clean workspace
rm(list = ls())


dat = read.table(
  file   = "Data/myDATA.csv",
  header = TRUE, 
  stringsAsFactors = FALSE,
  sep = ",",
  dec = ".",
  colClasses = rep('numeric',20),
  numerals = "no.loss"
)

indStart = which.max(dat[,1] >= 195001 )
indEnd   = nrow(dat)


Predictors_name = c("log_DP",    #   Log dividend-price ratio
                    "log_DY",    #   Log dividend yield
                    "log_EP",    #   Log earnings-price ratio
                    "log_DE",    #   Log dividend-payout ratio
                    "SVAR",      #   Stock variance
                    "BM",        #   Book-to-market ratio
                    "NTIS",      #   Net equity expansion
                    "TBL",       #   Treasury bill rate
                    "LTY",       #   Long-term yield
                    "LTR",       #   Long-term return
                    "TMS",       #   Term spread
                    "DFY",       #   Default yield spread
                    "DFR",       #   Default return spread
                    "INFL",      #   Inflation
                    "RET",       #   Last return
                    "RET2"       #   Last return squared
)

# Allocate memory for predictors 
Predictors = matrix( numeric(), 
                     nrow = nrow(dat), 
                     ncol = length(Predictors_name), 
                     dimnames = list(c(dat[,1]),Predictors_name ))

# Log dividend-price ratio and log dividend yield
Predictors[,"log_DP"] = log(dat$D12/dat$Index)
Predictors[13:nrow(dat),"log_DY"] = log(dat$D12[13:nrow(dat)]/dat$Index[1:(nrow(dat)-12)])

# Log earnings-price ratio and log dividend-payout ratio
Predictors[,"log_EP"] = log(dat$E12 / dat$Index)
Predictors[,"log_DE"] = log(dat$D12 / dat$E12)

# Predictors that do not require any transformation
Predictors[,"SVAR"]   = dat$svar
Predictors[,"BM"]     = dat$b.m
Predictors[,"NTIS"]   = dat$ntis
Predictors[,"TBL"]    = dat$tbl
Predictors[,"LTY"]    = dat$lty
Predictors[,"LTR"]    = dat$ltr

# Term spread
Predictors[, "TMS"]   = dat$lty - dat$tbl

# Default Yield spread
Predictors[,"DFY"]    = dat$BAA   - dat$AAA

# Default Return spread
Predictors[,"DFR"]    = dat$corpr - dat$ltr

# Inflation (delayed by one month to account for delay in CPI release)
Predictors[2:nrow(dat),"INFL"]   = dat$infl[-nrow(dat)]

# Lagged Return
Predictors[2:nrow(dat),"RET"]   = dat$infl[-nrow(dat)]

# Lagged Return Squared
Predictors[2:nrow(dat),"RET2"]   = dat$infl[-nrow(dat)]^2


# Compute Market Excess Return
ER = dat$Index.total.return - dat$Rfree

# Parameters
riskAversion = c(3, 5)
alpha        = c(0.01, 0.50, 0.85)
PortfolioReturn      = list()



for(ii in Predictors_name){
  PortfolioReturn[[paste("LM",ii, sep = "-")]] =   matrix( numeric(), 
                                                           nrow = nrow(dat), 
                                                           ncol = length(riskAversion), 
                                                           dimnames = list( c(dat[,1]),
                                                                            paste("RA",riskAversion, sep = "-") ))
}

PortfolioReturn[[paste("LM","ALL", sep = "-")]] =   matrix( numeric(), 
                                                            nrow = nrow(dat), 
                                                            ncol = length(riskAversion), 
                                                            dimnames = list( c(dat[,1]),
                                                                             paste("RA",riskAversion, sep = "-") ))
PortfolioReturn[[paste("LASSO","ALL", sep = "-")]] =   matrix( numeric(), 
                                                               nrow = nrow(dat), 
                                                               ncol = length(riskAversion), 
                                                               dimnames = list( c(dat[,1]),
                                                                                paste("RA",riskAversion, sep = "-") ))

PortfolioReturn[[paste("RIDGE","ALL", sep = "-")]] =   matrix( numeric(), 
                                                               nrow = nrow(dat), 
                                                               ncol = length(riskAversion), 
                                                               dimnames = list( c(dat[,1]),
                                                                                paste("RA",riskAversion, sep = "-") ))

PortfolioReturn[[paste("A.LASSO","ALL", sep = "-")]] =   matrix( numeric(), 
                                                                 nrow = nrow(dat), 
                                                                 ncol = length(riskAversion), 
                                                                 dimnames = list( c(dat[,1]),
                                                                                  paste("RA",riskAversion, sep = "-") ))


for(A in alpha ){
  PortfolioReturn[[paste0("ENET_",round(A,digits = 2),"-ALL", sep = "-")]] =   matrix( numeric(),
                                                                                       nrow = nrow(dat),
                                                                                       ncol = length(riskAversion),
                                                                                       dimnames = list( c(dat[,1]),
                                                                                                        paste("RA",riskAversion, sep = "-") ))  
}

###########
MLucas =   matrix( numeric(), 
                                                                 nrow = length(indStart:indEnd), 
                                                                 ncol = 1, 
                                                                 dimnames = list( c(dat[indStart:indEnd,1]),
                                                                                   'RSE'))
##########

# For loop of all the months from January 1950 to December 2015
# Progress bar
pb   <- txtProgressBar(indStart, indEnd, style=3)
TIME <- Sys.time()
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
for(t in indStart:indEnd){
  setTxtProgressBar(pb, t)
  # Dates available to estimate risk premium and variance
  #indDates = 2:(t-1)
  indDates = (t-501):(t-1)
  # Compute the sample variance on the sample of
  # returns available (up to month t - 1)
  # sampleVar = var(ER[indDates], na.rm = TRUE)
  sampleVar = var(dat$Index.total.return[indDates], na.rm = TRUE)
  #############################################
  # OLS: Naive Linear Model
  #############################################
  
  # Considering just one Predictor
  for(pred_name in Predictors_name){
    # Find the sample for which both returns and the 
    # predictor (lagged one period) is available
    AvailableData =  na.omit(
      cbind(Predictors[indDates-1,pred_name], ER[indDates] )
    )
    colnames(AvailableData) <- c(pred_name,"ER")
    AvailableData = data.frame(AvailableData)
    # Estimate the linear regresssion with the lagged predictor
    
    # Risk Premium Forecast
    LM = lm( as.formula(as.formula(paste0("ER ~ ",paste(pred_name, collapse=" + ")))) , 
             data = AvailableData)
    df.pred = as.data.frame( t(Predictors[t-1,pred_name]) ) 
    colnames(df.pred) <- pred_name
    RP_Forecast = predict(LM,  df.pred   )
    for(RA in riskAversion){
      PortfolioReturn[[paste("LM",pred_name, sep = "-")]][t, paste("RA",RA, sep = "-")] = 
        dat$Rfree[t]  + RP_Forecast/(RA*sampleVar) * ER[t]
    } 
  } 
  
  # Considering ALL predictors
  
  # Find the sample for which both returns and the 
  # predictor (lagged one period) is available
  pred_name = Predictors_name
  
  AvailableData =  na.omit(
    cbind(Predictors[indDates-1,pred_name], ER[indDates] )
  )
  colnames(AvailableData) <- c(pred_name,"ER")
  AvailableData = data.frame(AvailableData)
  # Estimate the linear regresssion with the lagged predictor
  
  # Risk Premium Forecast
  LM = lm( as.formula(as.formula(paste0("ER ~ ",paste(pred_name, collapse=" + ")))) , 
           data = AvailableData)
  df.pred = as.data.frame( t(Predictors[t-1,pred_name]) ) 
  colnames(df.pred) <- pred_name
  RP_Forecast = predict(LM,  df.pred  )
  
  for(RA in riskAversion){
    PortfolioReturn[[paste("LM","ALL", sep = "-")]][t, paste("RA",RA, sep = "-")] = 
      dat$Rfree[t]  + RP_Forecast/(RA*sampleVar) * ER[t]
  }
  
  
  #############################################
  # LASSO
  #############################################
  
  # Considering ALL predictorsr
  pred_name = Predictors_name
  AvailableData =  na.omit(
    cbind(Predictors[indDates-1,pred_name], ER[indDates] )
  )
  colnames(AvailableData) <- c(pred_name,"ER")
  
  # Estimate the LASSO regresssion with the lagged predictor
  cv.LASSO       = cv.glmnet( AvailableData[,pred_name], 
                              AvailableData[,"ER"],
                              parallel = TRUE)
  
  # Risk Premium Forecast
  RP_Forecast = predict(cv.LASSO,  newx = t(Predictors[t-1,pred_name]), s = "lambda.1se"  )
  
  # Portfolio Performance
  for(RA in riskAversion){
    PortfolioReturn[[paste("LASSO","ALL", sep = "-")]][t, paste("RA",RA, sep = "-")] = 
      dat$Rfree[t]  + RP_Forecast/(RA*sampleVar) * ER[t]
  }
  
  #############################################
  # RIDGE
  #############################################
  
  # Considering ALL predictorsr
  pred_name = Predictors_name
  AvailableData =  na.omit(
    cbind(Predictors[indDates-1,pred_name], ER[indDates] )
  )
  colnames(AvailableData) <- c(pred_name,"ER")
  
  # Estimate the LASSO regresssion with the lagged predictor
  cv.RIDGE       = cv.glmnet( AvailableData[,pred_name], 
                              AvailableData[,"ER"],
                              alpha    = 0,            # RIDGE
                              parallel = TRUE)  
  
  # Risk Premium Forecast
  RP_Forecast = predict(cv.RIDGE,  newx = t(Predictors[t-1,pred_name]), s = "lambda.1se"  )
  
  # Portfolio Performance
  for(RA in riskAversion){
    PortfolioReturn[[paste("RIDGE","ALL", sep = "-")]][t, paste("RA",RA, sep = "-")] = 
      dat$Rfree[t]  + RP_Forecast/(RA*sampleVar) * ER[t]
  }
  
  
  #############################################
  # A.LASSO
  #############################################
  gamma = 1
  W_A.LASSO = 1/abs(matrix( coef(cv.RIDGE, s=cv.RIDGE$lambda.min)     ) )^gamma 
  # Considering ALL predictorsr
  pred_name = Predictors_name
  AvailableData =  na.omit(
    cbind(Predictors[indDates-1,pred_name], ER[indDates] )
  )
  colnames(AvailableData) <- c(pred_name,"ER")
  
  # Estimate the LASSO regresssion with the lagged predictor
  cv.A.LASSO <- cv.glmnet( AvailableData[,pred_name], 
                           AvailableData[,"ER"], 
                           alpha          = 1, 
                           penalty.factor = W_A.LASSO,
                           parallel       = TRUE)
  
  # Risk Premium Forecast
  RP_Forecast = predict(cv.A.LASSO,  newx = t(Predictors[t-1,pred_name]), s = "lambda.1se"  )
  
  # Portfolio Performance
  for(RA in riskAversion){
    PortfolioReturn[[paste("A.LASSO","ALL", sep = "-")]][t, paste("RA",RA, sep = "-")] = 
      dat$Rfree[t]  + RP_Forecast/(RA*sampleVar) * ER[t]
  }
  
  for(A in alpha){
    #############################################
    # ELASTIC NET
    #############################################
    
    # Considering ALL predictorsr
    pred_name = Predictors_name
    AvailableData =  na.omit(
      cbind(Predictors[indDates-1,pred_name], ER[indDates] )
    )
    colnames(AvailableData) <- c(pred_name,"ER")
    
    # Estimate the LASSO regresssion with the lagged predictor
    cv.ENET       = cv.glmnet(  AvailableData[,pred_name], 
                                AvailableData[,"ER"],
                                alpha    = A,   
                                parallel = TRUE)  
    
    # Risk Premium Forecast
    RP_Forecast = predict(cv.ENET,  newx = t(Predictors[t-1,pred_name]), s = "lambda.1se"  )
    
    # Portfolio Performance
    for(RA in riskAversion){
      PortfolioReturn[[paste0("ENET_",round(A,digits = 2),"-ALL", sep = "-")]][t, paste("RA",RA, sep = "-")] = 
        dat$Rfree[t]  + RP_Forecast/(RA*sampleVar) * ER[t]
    }
    
  }
  
  
}


# options(warn = oldw)
stopCluster(cl)

# Compute realized certainty equivalents



CE =  matrix( numeric(), 
              nrow = length(PortfolioReturn), 
              ncol = length(riskAversion), 
              dimnames = list( names(PortfolioReturn),
                               paste("RA",riskAversion, sep = "-") ))

for(i in names(PortfolioReturn)){
  for(RA in riskAversion){
    R_aux = PortfolioReturn[[i]][indStart:indEnd,paste("RA",RA, sep = "-") ]
    CE[i,paste("RA",RA, sep = "-") ] = mean(R_aux ,na.rm = TRUE) - 0.5*RA*var(R_aux, na.rm = TRUE)
  }
}


print(CE*12*100, digits = 3)
print(Sys.time() - TIME)

# my.dat = cbind(Predictors,ER)
# my.dat = na.omit(my.dat)
# 
# ER = my.dat[,"ER"]
# Predictors = my.dat[, Predictors_name]
# 
# fit = glmnet(Predictors, ER)
# 
# plot(fit, label = TRUE)
# 
# plot(fit, xvar = "norm",   label = TRUE)
# plot(fit, xvar = "lambda", label = TRUE)
# plot(fit, xvar = "dev",    label = TRUE)
# 
# print(fit)
# cv.fit = cv.glmnet(Predictors, ER, type.measure = "mae", alpha = 0.15)
# plot(cv.fit)
# cv.fit$lambda.1se


