library(forecast)
library(xgboost)
library(tidyverse)
library(Mcomp)
library(bsts)
library(parallel)
library(doMC)

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
  model <- bsts(ts, state.specification = ss, niter = 2500, family = "gaussian")
  bsts <- predict(model, horizon = h, burn = 200)$mean
  
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

medianEnsemble <- function(ts, h){
  # ts is a time series object
  # h is the forecast horizon, the number of data points into the future to forecast
  
  # Fit the automatic models to the ts object
  # Arima Forecast
  aa <- forecast(auto.arima(ts), h = h)$mean
  
  # Neural Network
  ss <- AddLocalLinearTrend(list(), ts)
  model <- bsts(ts,state.specification = ss, niter = 1000)
  bsts <- predict(model, horizon = h)$mean
  
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


xgStackEnsemble <- function(ts, h, validationSize = round(length(ts)/2)){
  # Demonstration of Ensemble Stacking Technique for Time series forecasts
  
  # Divide data into train and validation sets
  train <- head(ts, length(ts)-validationSize)
  valid <- tail(ts, validationSize)
  
  # Using training data, make prediction on validation set using nEnsemble
  valMat <- naiveEnsemble(train, validationSize)[[1]]
  cat("Validation Predications Made")
  
  # Use XGBoost algorithm to find best nonlinear transformation of forecasted validation
  # data points into actual validation observations
  xgStack <- xgboost(data = as.matrix(valMat), 
                     label = as.numeric(valid),
                     nrounds = 2000,
                     eta = 0.01)
  cat("XGB Stack completed")
  
  # Combine train and validation sets to make forecasts h steps into the horizon
  fcMat <- naiveEnsemble(ts, h)[[1]]
  cat("Forecast Matrix Created")
  # Apply the learned XGBoost nonlinear transformation
  xgPred <- predict(xgStack, fcMat)
  cat("Predictions Made")
  # Return Predictions
  return(as.numeric(xgPred))
}

linStackEnsemble <- function(ts, h, validationSize = round(length(ts)/2)){
  # Demonstration of Ensemble Stacking Technique for Time series forecasts
  
  # Divide data into train and validation sets
  train <- head(ts, length(ts)-validationSize)
  valid <- tail(ts, validationSize)
  
  # Using training data, make prediction on validation set using nEnsemble
  valMat <- naiveEnsemble(train, validationSize)[[1]]
  cat("Validation Predications Made")
  
  valDF <- as.data.frame(valMat)
  names(valDF) <- c("aa", "bsts", "ets", "theta")
  valDF$y <- as.numeric(valid)
  
  # Use XGBoost algorithm to find best nonlinear transformation of forecasted validation
  # data points into actual validation observations
  linStack <- loess(y ~ aa + bsts + ets + theta,
                 data = valDF,
                 control = loess.control(surface = "direct"))
  cat("Linear Stack completed")
  
  # Combine train and validation sets to make forecasts h steps into the horizon
  fcMat <- naiveEnsemble(ts, h)[[1]]
  fcDF <- as.data.frame(fcMat)
  names(fcDF) <- c("aa", "bsts", "ets", "theta")
  cat("Forecast Matrix Created")
  
  # Apply the linear transformation
  linPred <- predict(linStack, fcDF)
  cat("Predictions Made")
  # Return Predictions
  return(as.numeric(linPred))
}


testing <- function(M3obj){
  mean <- mean(c(M3obj$x, M3obj$xx))
  ts <- M3obj$x
  xx <- M3obj$xx
  h <- length(xx)
  Period <- M3obj$eriod
  Series <- M3obj$sn
  seasonal <- findfrequency(ts)
  
  aa <- forecast(auto.arima(ts), h = h)$mean
  
  ss <- AddLocalLinearTrend(list(), ts)
  if(seasonal>1){ss <- AddSeasonal(ss, ts, nseasons = seasonal)}
  model <- bsts(ts, state.specification = ss, niter = 2500, family = "gaussian")
  bsts <- predict(model, horizon = h, burn = 200)$mean
  
  ets <- forecast(ets(ts), h = h)$mean
  
  theta <- forecast(thetaf(ts, h = h), h = h)$mean
  
  fcMat <- as.matrix(cbind(as.numeric(aa), as.numeric(bsts), as.numeric(ets), 
                           as.numeric(theta)))
  
  BEAT <-rowMeans(fcMat)
  
  mBEAT <- apply(fcMat, 1, median)
  
  XGStack <- xgStackEnsemble(ts, h = h)
  
  
  fcList <- list(accuracy(aa,xx), accuracy(bsts,xx), accuracy(ets,xx), 
                 accuracy(theta,xx), accuracy(BEAT,xx), accuracy(mBEAT,xx),
                 accuracy(XGStack, xx))
  
  out <- as.data.frame(do.call(rbind,fcList))
  out$Series <- Series
  out$Period <- Period
  out$Method <- c("Auto.Arima", "BSTS", "ETS", "THETA", "BEAT", "mBEAT", "XGStack")
  rownames(out) <- NULL
  out$ME <- out$ME/mean
  out$RMSE <- out$RMSE/mean
  out$MAE <- out$MAE/mean
  out %>% dplyr::select(Series, Period, Method, ME, RMSE, MAE, MAPE)
}


# # Non Parallel
# list <- list()
# for(i in 1:100){
#   cat(i, '\n')
#   out <- testing(M3[[i]])
#   list[[i]] <- out
# }
# as.data.frame(do.call(rbind,list)) %>% group_by(Method) %>% summarise(RMSE = mean(RMSE), ME = mean(ME), MAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy
# 
# list <- list()
# n = 1
# foreach(i=1:length(M3)) %dopar% {
#   cat(n, "\n")
#   out <- testing(M3[[i]])
#   list[[i]] <- out
#   n = n + 1
# }



#samp <- sample(1:3003, 100)
timeOut <- system.time({ 
  outDF <- foreach(i = 1:2000) %dopar% {
    out <- testing(M3[[i]])
  }
})

chunk <- as.data.frame(do.call(rbind, outDF))
final <- chunk
final <- rbind(final, chunk)

# run at end
final %>% group_by(Method) %>% summarise(RMSE = mean(RMSE), ME = mean(ME), MAE = mean(MAE), MAPE = mean(MAPE)) -> accuracy

