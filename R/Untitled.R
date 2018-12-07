library(forecast)
library(randomForest)
library(xgboost)

nEnsemble <- function(ts, horizon){
  
  # Arima Forecast
  aa <- forecast(auto.arima(ts, approximation = F, stepwise = F), h = horizon)$mean
  
  # Neural Network
  nn <- forecast(nnetar(ts, scale.inputs = T, decay = .2), h = horizon)$mean
  
  # Exponential Smoothing
  ets <- forecast(ets(ts), h = horizon)$mean
  
  # Theta forecast
  theta <- forecast(thetaf(ts, h = horizon), h = horizon)$mean
  
  # Seasonal Naive
  snaive <- snaive(ts, h = horizon)$mean
  
  # Forecast Matrix
  fcMat <- as.matrix(cbind(as.numeric(aa), as.numeric(nn), as.numeric(ets), 
                           as.numeric(theta), as.numeric(snaive)))

  out <- list(predictionMatrix = fcMat, forecast = rowMeans(fcMat))

  return(out)
}

rfStackEnsemble <- function(ts, horizon, validationSize = 30){
  
  train <- head(ts, length(ts)-validationSize)
  valid <- tail(ts, validationSize)
  
  # Find Validatin Predictions
  
  valMat <- nEnsemble(train, validationSize)[[1]]
  cat("Validation Predications Made")
  rfStack <- randomForest(x = valMat, y = as.numeric(valid), ntree = 100)
  cat("RF Stack completed")
  fcMat <- nEnsemble(ts, horizon)[[1]]
  cat("FC Matrix Created")
  rfPred <- predict(rfStack, fcMat)
  cat("Predictions Made")
  return(rfPred)
}

xgStackEnsemble <- function(ts, horizon, validationSize = 30){
  
  train <- head(ts, length(ts)-validationSize)
  valid <- tail(ts, validationSize)
  
  # Find Validatin Predictions
  
  valMat <- nEnsemble(train, validationSize)[[1]]
  print(valMat)
  print(as.numeric(valid))
  cat("Validation Predications Made")
  xgStack <- xgboost(data = as.matrix(valMat), 
                     label = as.numeric(valid),
                     nrounds = 200,
                     eta = 0.1)
  cat("RF Stack completed")
  fcMat <- nEnsemble(ts, horizon)[[1]]
  cat("FC Matrix Created")
  xgPred <- predict(xgStack, fcMat)
  cat("Predictions Made")
  return(xgPred)
}
