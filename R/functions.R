library(forecast)
library(xgboost)
library(nnet)

Ensemble <- function(ts, horizon){
  # Auto Arima
  aa <- forecast(auto.arima(ts, approximation = F, stepwise = F), h = horizon)$mean
  
  # Neural Net
  nn <- forecast(nnetar(ts, scale.inputs = T, decay = .5), h = horizon)$mean
  
  # Exponential Smoothing
  ets <- forecast(ets(ts), h = horizon)$mean
  
  # Theta Method
  theta <- forecast(thetaf(ts, h = horizon), h = horizon)$mean
  
  # Seasonal Naive Forecast
  snaive <- snaive(ts, h = horizon)$mean
  
  fcMat <- as.matrix(cbind(aa, nn, ets, theta, snaive))
  
  return(fcMat)
}

nEnsemble <- function(ts, horizon){
  fcMat <- Ensemble(ts, horizon)
  forecast <- as.numeric(rowMeans(fcMat))
  return(forecast)
}

xgEnsemble <- function(ts, horizon, validationSize, nrounds = 500, params = list(eta = 0.001)){
  train <- head(ts, length(ts) - validationSize)
  validation <- tail(ts, validationSize)
  
  validationMatrix <- Ensemble(train, validationSize)
  xgStack <- xgboost(data = validationMatrix,
                     label = as.numeric(validation),
                     nrounds = nrounds,
                     params = params)
  
  predictionMat <- Ensemble(ts, horizon)
  forecast <- predict(xgStack, predictionMat)
  return(forecast)
}

ensembleAccuracy <- function(fcMat, testSet){
  return(apply(fcMat, 2, function(x) accuracy(x, testSet)))
}
