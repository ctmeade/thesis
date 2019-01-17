library(forecast)
library(xgboost)
library(tidyverse)
library(Mcomp)
library(bsts)
library(SuperLearner)
library(parallel)
library(doMC)
library(e1071)
library(bartMachine)

doMC::registerDoMC(cores = detectCores())

naiveEnsemble <- function(ts, h){
  # ts is a time series object
  # h is the forecast horizon, the number of data points into the future to forecast
  
  # Fit the automatic models to the ts object
  # Arima Forecast
  aa <- forecast(auto.arima(ts), h = h)$mean
  
  # Neural Network
  seasonal <- findfrequency(ts)
  ss <- AddLocalLinearTrend(list(), ts)
  if(seasonal>1){ss <- AddSeasonal(ss, ts, nseasons = seasonal)}
  model <- bsts(ts, state.specification = ss, niter = 1000, family = "gaussian")
  bsts <- predict(model, horizon = h, burn = 100)$mean
  
  # Exponential Smoothing
  ets <- forecast(ets(ts), h = h)$mean
  
  # Theta forecast
  theta <- forecast(thetaf(ts, h = h), h = h)$mean
  
  # Combine each h by 1 forecast into a h by p matrix, where p is the number
  # of models fit to the data
  
  # Forecast Matrix
  fcMat <- as.matrix(cbind(as.numeric(aa), as.numeric(bsts), as.numeric(ets), 
                           as.numeric(theta)))
  
  # Return the Forecast Matrix and the average forecast at each point on
  # the horizon as a list object
  out <- list(predictionMatrix = fcMat, forecast = rowMeans(fcMat))

  return(out)
}

EATEnsemble <- function(ts, h){
  # ts is a time series object
  # h is the forecast horizon, the number of data points into the future to forecast
  
  # Fit the automatic models to the ts object
  # Arima Forecast
  aa <- forecast(auto.arima(ts), h = h)$mean
  
  # Exponential Smoothing
  ets <- forecast(ets(ts), h = h)$mean
  
  # Theta forecast
  theta <- forecast(thetaf(ts, h = h), h = h)$mean
  
  # Combine each h by 1 forecast into a h by p matrix, where p is the number
  # of models fit to the data
  
  # Forecast Matrix
  fcMat <- as.matrix(cbind(as.numeric(aa), as.numeric(ets), 
                           as.numeric(theta)))
  
  # Return the Forecast Matrix and the average forecast at each point on
  # the horizon as a list object
  out <- list(predictionMatrix = fcMat, forecast = rowMeans(fcMat))
  
  return(out)
}

medianEnsemble <- function(ts, h){
  # ts is a time series object
  # h is the forecast horizon, the number of data points into the future to forecast
  
  # Fit the automatic models to the ts object
  # Arima Forecast
  aa <- forecast(auto.arima(ts), h = h)$mean
  
  # Neural Network
  seasonal <- findfrequency(ts)
  ss <- AddLocalLinearTrend(list(), ts)
  if(seasonal>1){ss <- AddSeasonal(ss, ts, nseasons = seasonal)}
  model <- bsts(ts, state.specification = ss, niter = 1000, family = "gaussian")
  bsts <- predict(model, horizon = h, burn = 100)$mean
  
  # Exponential Smoothing
  ets <- forecast(ets(ts), h = h)$mean
  
  # Theta forecast
  theta <- forecast(thetaf(ts, h = h), h = h)$mean
  
  # Combine each h by 1 forecast into a h by p matrix, where p is the number
  # of models fit to the data
  
  # Forecast Matrix
  fcMat <- as.matrix(cbind(as.numeric(aa), as.numeric(bsts), as.numeric(ets), 
                           as.numeric(theta)))
  
  # Return the Forecast Matrix and the average forecast at each point on
  # the horizon as a list object
  out <- list(predictionMatrix = fcMat, forecast = apply(fcMat, 1, median))
  
}

nonparallel_baggedBEAT <- function(ts, h){
  bootList <- bld.mbb.bootstrap(ts, 10)
  
  out <- list()
  for(i in 1:length(bootList)){
    cat(i)
    out[[i]] <- as.numeric(naiveEnsemble(bootList[[i]], h = h)[[2]])
  }
  as.data.frame(do.call(rbind, out))
}

error.resamp <- function(x, num, block_size=NULL) {
  freq <- frequency(x)
  if (is.null(block_size)) {
    block_size <- ifelse(freq > 1, 2 * freq, min(8, floor(length(x) / 2)))
  }
  
  xs <- list()
  xs[[1]] <- x # the first series is the original one
  
  if (num > 1) {
    # Box-Cox transformation
    if (min(x) > 1e-6) {
      lambda <- BoxCox.lambda(x, lower = 0, upper = 1)
    } else {
      lambda <- 1
    }
    x.bc <- BoxCox(x, lambda)
    lambda <- attr(x.bc, "lambda")
    
    if (freq > 1) {
      # STL decomposition
      x.stl <- stl(ts(x.bc, frequency = freq), "per")$time.series
      seasonal <- x.stl[, 1]
      trend <- x.stl[, 2]
      remainder <- x.stl[, 3]
    } else {
      # Loess
      trend <- 1:length(x)
      suppressWarnings(x.loess <- loess(x.bc ~ trend, span = 6 / length(x), degree = 1))
      seasonal <- rep(0, length(x))
      trend <- x.loess$fitted
      remainder <- x.loess$residuals
    }
    
    # Bootstrap some series, using MBB
    for (i in 2:num) {
      xs[[i]] <- InvBoxCox(trend + seasonal + rnorm(length(remainder), 0, sd(remainder)), lambda)
    }
  }
  
  xs
}

nonparallel_pBEAT <- function(ts, h){
  bootList <- bld.mbb.bootstrap(ts, 10)
  
  out <- list()
  for(i in 1:length(bootList)){
    cat(i)
    out[[i]] <- as.numeric(naiveEnsemble(bootList[[i]], h = h)[[2]])
  }
  as.data.frame(do.call(rbind, out))
}

baggedBEAT <- function(ts, h){
  bootList <- bld.mbb.bootstrap(ts, 10)
  
  outDF <- foreach(i = 1:length(bootList)) %dopar% {
    out <- as.numeric(naiveEnsemble(bootList[[i]], h = h)[[2]])
  }
  
  as.data.frame(do.call(rbind, outDF))
  
}

baggedEAT <- function(ts, h){
  bootList <- bld.mbb.bootstrap(ts, 10)
  
  outDF <- foreach(i = 1:length(bootList)) %dopar% {
    out <- as.numeric(EATEnsemble(bootList[[i]], h = h)[[2]])
  }
  
  as.data.frame(do.call(rbind, outDF))
  
}

pBEAT <- function(ts, h){
  bootList <- error.resamp(ts, 10)
  
  outDF <- foreach(i = 1:length(bootList)) %dopar% {
    out <- as.numeric(naiveEnsemble(bootList[[i]], h = h)[[2]])
  }
  
  as.data.frame(do.call(rbind, outDF))
}

testing <- function(M3obj){
  mean <- mean(M3obj$x)
  ts <- M3obj$x
  xx <- M3obj$xx
  h <- length(xx)
  Period <- M3obj$period
  Series <- M3obj$sn
  seasonal <- frequency(ts)
  
  aa <- forecast(auto.arima(ts), h = h)$mean
  
  ss <- AddLocalLinearTrend(list(), ts)
  if(seasonal>1){ss <- AddSeasonal(ss, ts, nseasons = seasonal)}
  model <- bsts(ts, state.specification = ss, niter = 1000, family = "gaussian")
  bsts <- predict(model, horizon = h, burn = 100)$mean
  
  ets <- forecast(ets(ts), h = h)$mean
  
  theta <- forecast(thetaf(ts, h = h), h = h)$mean
  
  fcMat <- as.matrix(cbind(as.numeric(bsts), as.numeric(ets), as.numeric(aa), 
                           as.numeric(theta)))
  
  BEA <- rowMeans(fcMat[,c(1,2,3)])
  EAT <- rowMeans(fcMat[,c(2,3,4)])
  BAT <- rowMeans(fcMat[,c(1,3,4)])
  BET <- rowMeans(fcMat[,c(1,2,4)])
  
  BEAT <-rowMeans(fcMat)
  
  mBEAT <- apply(fcMat, 1, median)
  
  baggedBEAT <- baggedBEAT(ts, h)
  meanBaggedBEAT <- as.numeric(colMeans(baggedBEAT))
  medianBaggedBEAT <- apply(baggedBEAT, 2, median)

  pBEAT <- pBEAT(ts, h)
  meanPertBEAT <- as.numeric(colMeans(pBEAT))
  medianPertBEAT <- apply(pBEAT, 2, median)
  
  fcList <- list(accuracy(aa,xx), 
                 accuracy(bsts,xx), 
                 accuracy(ets,xx), 
                 accuracy(theta,xx), 
                 accuracy(BEA,xx), 
                 accuracy(EAT,xx), 
                 accuracy(BAT,xx), 
                 accuracy(BET,xx), 
                 accuracy(BEAT,xx), 
                 accuracy(mBEAT,xx),
                 accuracy(meanBaggedBEAT, xx),
                 accuracy(medianBaggedBEAT, xx),
                 accuracy(meanPertBEAT, xx),
                 accuracy(medianPertBEAT, xx)
  )
  
  out <- as.data.frame(do.call(rbind,fcList))
  out$Series <- Series
  out$Period <- Period
  out$Method <- c("Auto.Arima", "BSTS", "ETS", "THETA", "BEA", "EAT", "BAT", "BET", "BEAT", "medianBEAT", "meanBaggedBEAT", "medianBaggedBEAT", "meanPertBEAT", "medianPertBEAT")
  rownames(out) <- NULL
  out$ME <- out$ME/mean
  out$RMSE <- out$RMSE/mean
  out$MAE <- out$MAE/mean
  out %>% dplyr::select(Series, Period, Method, RMSE, MAE, MAPE)
}

testing2 <- function(M3obj){
  mean <- mean(M3obj$x)
  ts <- M3obj$x
  xx <- M3obj$xx
  h <- length(xx)
  Period <- M3obj$period
  Series <- M3obj$sn
  
  baggedEAT <- baggedEAT(ts, h)
  meanBaggedEAT <- as.numeric(colMeans(baggedEAT))
  
  fcList <- list(accuracy(meanBaggedEAT,xx)
                 
  )
  
  out <- as.data.frame(do.call(rbind,fcList))
  out$Series <- Series
  out$Period <- Period
  out$Method <- c("meanBaggedEAT")
  rownames(out) <- NULL
  out$ME <- out$ME/mean
  out$RMSE <- out$RMSE/mean
  out$MAE <- out$MAE/mean
  out %>% dplyr::select(Series, Period, Method, RMSE, MAE, MAPE)
}

samp <- sample(1:3003, 1)
timeOut <- system.time({ 
  outDF <- foreach(i = 1:3003) %dopar% {
    out <- testing2(M3[[i]])
  }
})

chunk <- as.data.frame(do.call(rbind, outDF))
final <- chunk
final <- rbind(final, chunk)
final %>% group_by(Period, Method) %>% summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy



seqList <- split(1:3003, ceiling(seq_along(1:3003)/10))
for(i in 1:length(seqList)){
  outDF <- foreach(i = seqList[[i]]) %dopar% {
    out <- testing2(M3[[i]])
  }
  chunk <- as.data.frame(do.call(rbind, outDF))
  write.csv(chunk, paste(i, ".csv", sep = ""))
}


filenames <- list.files(full.names=TRUE)
All <- lapply(filenames,function(i){
  read.csv(i, header=T, skip=0)
})
data <- as.data.frame(do.call(rbind, All))

data %>% 
  group_by(Method) %>% 
  summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy

data %>% 
  dplyr::filter(Period == "MONTHLY") %>% 
  group_by(Method) %>% 
  summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> monthly

data %>% 
  dplyr::filter(Period == "YEARLY") %>% 
  group_by(Method) %>% 
  summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> yearly

data %>% 
  dplyr::filter(Period == "QUARTERLY") %>% 
  group_by(Method) %>% 
  summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> quarterly

data %>% 
  dplyr::filter(Period == "OTHER") %>% 
  group_by(Method) %>% 
  summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> other
# # run at end
# final %>% group_by(Period, Method) %>% summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy
# 
# M3df <- list()
# for(i in 1:length(M3)){
#   M3df[[i]] <- list(Name = M3[[i]]$sn, Period = M3[[i]]$period)
# }
# M3df <- as.data.frame(do.call(rbind,M3df))
# M3df <- data.frame(as.character(M3df$Name), as.character(M3df$Period))
# names(M3df) <- c("Series", "Period")                      
# 
# final %>% left_join(M3df, by = c("Series")) %>% group_by(Period, Method) %>% summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy
# 
# testing2 <- function(M3obj){
#   mean <- mean(c(M3obj$x, M3obj$xx))
#   ts <- M3obj$x
#   xx <- M3obj$xx
#   h <- length(xx)
#   Period <- M3obj$period
#   Series <- M3obj$sn
#   bBEAT <- baggedBEAT(ts, h)
#   
#   fcList <- list(accuracy(bBEAT,h))
#   
#   out <- as.data.frame(do.call(rbind,fcList))
#   out$Series <- Series
#   out$Period <- Period
#   out$Method <- c("baggedBEAT")
#   rownames(out) <- NULL
#   out$ME <- out$ME/mean
#   out$RMSE <- out$RMSE/mean
#   out$MAE <- out$MAE/mean
#   out %>% dplyr::select(Series, Period, Method, RMSE, MAE, MAPE)
# }
# 
# samp <- sample(1:3003, 100)
# timeOut <- system.time({ 
#   outDF <- foreach(i = samp) %dopar% {
#     out <- testing(M3[[i]])
#   }
# })
# 
# chunk <- as.data.frame(do.call(rbind, outDF))
# chunk %>% group_by(Method) %>% summarise(mRMSE = mean(RMSE), mMAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy
#
# xgStackEnsemble <- function(ts, h, validationSize = round(length(ts)/1.5)){
#   # Demonstration of Ensemble Stacking Technique for Time series forecasts
#   
#   # Divide data into train and validation sets
#   train <- head(ts, length(ts)-validationSize)
#   valid <- tail(ts, validationSize)
#   
#   # Using training data, make prediction on validation set using nEnsemble
#   valMat <- naiveEnsemble(train, validationSize)[[1]]
#   cat("Validation Predications Made")
#   
#   # Use XGBoost algorithm to find best nonlinear transformation of forecasted validation
#   # data points into actual validation observations
#   xgStack <- xgboost(data = as.matrix(valMat), 
#                      label = as.numeric(valid),
#                      nrounds = 3000,
#                      eta = 0.01)
#   cat("XGB Stack completed")
#   
#   # Combine train and validation sets to make forecasts h steps into the horizon
#   fcMat <- naiveEnsemble(ts, h)[[1]]
#   cat("Forecast Matrix Created")
#   # Apply the learned XGBoost nonlinear transformation
#   xgPred <- predict(xgStack, fcMat)
#   cat("Predictions Made")
#   # Return Predictions
#   return(as.numeric(xgPred))
# }
# superEnsemble <- function(ts, h, validationSize = round(length(ts)/1.5)){
#   # Demonstration of Ensemble Stacking Technique for Time series forecasts
#   
#   # Divide data into train and validation sets
#   train <- head(ts, length(ts)-validationSize)
#   valid <- tail(ts, validationSize)
#   
#   # Using training data, make prediction on validation set using nEnsemble
#   valMat <- as.data.frame(naiveEnsemble(train, validationSize)[[1]])
#   cat("Validation Predications Made")
#   
#   # Use XGBoost algorithm to find best nonlinear transformation of forecasted validation
#   # data points into actual validation observations
#   model <- SuperLearner(Y = as.numeric(valid), 
#                         X = valMat, 
#                         family = gaussian(),
#                         SL.library = c("SL.glmnet", "SL.xgboost", "SL.svm"))
#   cat("Super Stack completed")
#   
#   # Combine train and validation sets to make forecasts h steps into the horizon
#   fcMat <- naiveEnsemble(ts, h)[[1]]
#   fcMat <- as.data.frame(fcMat)
#   cat("Forecast Matrix Created")
#   # Apply the learned XGBoost nonlinear transformation
#   xgPred <- predict(model, fcMat)
#   cat("Predictions Made")
#   # Return Predictions
#   return(as.numeric(xgPred$pred))
# }
# linStackEnsemble <- function(ts, h, validationSize = round(length(ts)/2)){
#   # Demonstration of Ensemble Stacking Technique for Time series forecasts
#   
#   # Divide data into train and validation sets
#   train <- head(ts, length(ts)-validationSize)
#   valid <- tail(ts, validationSize)
#   
#   # Using training data, make prediction on validation set using nEnsemble
#   valMat <- naiveEnsemble(train, validationSize)[[1]]
#   cat("Validation Predications Made")
#   
#   valDF <- as.data.frame(valMat)
#   names(valDF) <- c("aa", "bsts", "ets", "theta")
#   valDF$y <- as.numeric(valid)
#   
#   # Use XGBoost algorithm to find best nonlinear transformation of forecasted validation
#   # data points into actual validation observations
#   linStack <- loess(y ~ aa + bsts + ets + theta,
#                  data = valDF,
#                  control = loess.control(surface = "direct"))
#   cat("Linear Stack completed")
#   
#   # Combine train and validation sets to make forecasts h steps into the horizon
#   fcMat <- naiveEnsemble(ts, h)[[1]]
#   fcDF <- as.data.frame(fcMat)
#   names(fcDF) <- c("aa", "bsts", "ets", "theta")
#   cat("Forecast Matrix Created")
#   
#   # Apply the linear transformation
#   linPred <- predict(linStack, fcDF)
#   cat("Predictions Made")
#   # Return Predictions
#   return(as.numeric(linPred))
# }


